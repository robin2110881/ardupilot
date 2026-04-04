#!/usr/bin/env python3
# This script was written mostly with the help of GitHub Copilot.
"""SymPy derivation helper for EKFGSF 5-state yaw model.

This is intentionally small and copy/paste friendly (unlike EKF3's full generator).
It derives:
- Discrete-time linearisation matrices F and G for x=[vN,vE,psi,pN,pE]
  with process-noise w=[dvx,dvy,dpsi] (front/right delta-velocity and yaw delta-angle).
- Optical-flow measurement model matching EKF3 optical flow fusion:
            z = ofDataDelayed.flowRadXYcomp  (rad/sec)
            losPred.x =  relVelSensor.y / range
            losPred.y = -relVelSensor.x / range
    where:
            relVelSensor = Tnb(roll,pitch,psi) * [vN, vE, vD] + omega_body x sensor_offset_body
    and roll/pitch, vD, omega_body, sensor offset, and range are treated as external inputs
    (consistent with GSF using external tilt/height/body-rate information).
- Position measurement model:
      z=[pN,pE],  h(x)=[pN,pE]

Run:
  python3 ekfgsf_yaw_5state_sympy.py > generated/ekfgsf_yaw_5state_autocode.txt

Notes:
- Optical flow sensors typically provide line-of-sight angular rates; in EKF3,
  gyro-compensation and frame transforms happen before fusion.
- If you want to account for height uncertainty (sigma_h^2) in the innovation variance,
    augment S by adding J*Cov_params*J^T where params can include height and tilt.
"""

from __future__ import annotations

import sympy as sp


def cheap_simplify(expr: sp.Expr) -> sp.Expr:
    """Fast-ish simplification suitable for codegen.

    Avoid SymPy's full `simplify()` which can explode on trig/rational expressions.
    """
    try:
        expr = sp.together(expr)
    except Exception:
        pass
    try:
        expr = sp.factor_terms(expr)
    except Exception:
        pass
    return expr


def simplify_linear_coeffs(expr: sp.Expr, *vars: sp.Symbol) -> sp.Expr:
    """Simplify coefficients of expressions linear in `vars`.

    If `expr` is affine (total degree <= 1) in the given vars, rewrite it as
    sum_i (ci * var_i) + c0 and simplify only the coefficients ci, c0.

    This is intended to keep optical-flow expressions tractable.
    """
    if not vars:
        return expr

    expr_in = expr
    try:
        expr_in = sp.together(expr_in)
    except Exception:
        pass

    try:
        poly = sp.Poly(expr_in, *vars, domain=sp.EX)
    except Exception:
        return expr

    try:
        if poly.total_degree() > 1:
            return expr
    except Exception:
        return expr

    out = sp.Integer(0)
    for monom, coef in poly.terms():
        coef_s = cheap_simplify(coef)
        term = coef_s
        for v, p in zip(vars, monom):
            if p == 0:
                continue
            if p == 1:
                term *= v
                continue
            # shouldn't happen for affine
            return expr
        out += term
    return out


def cexpr(expr: sp.Expr, *, simplifier=None) -> str:
    """C-ish expression string.

    IMPORTANT: Do not call full `simplify()` by default. Use `simplifier=` to
    optionally apply a cheap/local simplification pass.
    """
    if simplifier is not None:
        try:
            expr = simplifier(expr)
        except Exception:
            pass
    return sp.ccode(expr)


def emit_cse(
    assignments: list[tuple[str, sp.Expr]],
    temp_prefix: str = "t",
    *,
    simplifier=None,
) -> None:
    """Emit C++ assignments with common subexpression elimination.

    `assignments` is a list of (lhs, rhs).
    """
    rhs_list = [rhs for (_, rhs) in assignments]
    repl, reduced = sp.cse(rhs_list, symbols=sp.numbered_symbols(temp_prefix))
    for sym, expr in repl:
        print(f'const ftype {sym} = {cexpr(expr, simplifier=simplifier)};')
    for (lhs, _), rhs in zip(assignments, reduced):
        print(f'{lhs} = {cexpr(rhs, simplifier=simplifier)};')


def make_symm_P() -> tuple[sp.Matrix, list[list[sp.Symbol]]]:
    """Create a 5x5 symmetric covariance matrix with named elements P00..P44."""
    P = [[None for _ in range(5)] for _ in range(5)]
    syms: list[list[sp.Symbol]] = [[None for _ in range(5)] for _ in range(5)]  # type: ignore[assignment]
    for i in range(5):
        for j in range(i, 5):
            sym = sp.Symbol(f'P{i}{j}', real=True)
            P[i][j] = sym
            P[j][i] = sym
            syms[i][j] = sym
            syms[j][i] = sym
    return sp.Matrix(P), syms


def ekf_update_2d(x: sp.Matrix,
                 P: sp.Matrix,
                 h: sp.Matrix,
                 z: sp.Matrix,
                 R: sp.Symbol,
                 extra_S: sp.Matrix | None = None,
                 ) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """2D EKF update (pred - meas innovation), with isotropic R.

    Returns (innov(2x1), S(2x2), K(5x2), xnew(5x1), Pnew(5x5)).
    """
    H = h.jacobian(x)  # (2x5)

    innov = h - z  # EKF3 uses predicted - measured

    # PH^T and S
    PHt = P * H.T  # (5x2)

    #S = sp.simplify(H * PHt + R * sp.eye(2))  # (2x2)
    S = H * PHt + R * sp.eye(2)  # (2x2)
    if extra_S is not None:
        #S = sp.simplify(S + extra_S)
        S = S + extra_S

    # IMPORTANT: do NOT call S.inv() here.
    # SymPy's general symbolic inverse can be extremely slow even for 2x2 when the
    # entries are large (optical-flow model creates big trig/rational expressions).
    # Use explicit 2x2 adjugate/det inverse instead.
    detS = S[0, 0] * S[1, 1] - S[0, 1] * S[1, 0]
    S_inv = sp.Matrix(
        [
            [S[1, 1], -S[0, 1]],
            [-S[1, 0], S[0, 0]],
        ]
    ) / detS
    # K
    #K = sp.simplify(PHt * S_inv)
    K = PHt * S_inv

    # State update
    xnew = x - K * innov

    Pnew = P - PHt * S_inv * PHt.T
    return innov, S, K, xnew, Pnew


def main() -> None:
    # =========================
    # Discrete process model linearisation
    # =========================

    # state
    vN, vE, psi, pN, pE = sp.symbols('vN vE psi pN pE', real=True)

    # inputs/noise
    dvx, dvy, dpsi = sp.symbols('dvx dvy dpsi', real=True)
    dt = sp.symbols('dt', positive=True, real=True)

    s = sp.sin(psi)
    c = sp.cos(psi)

    # rotate front/right delta-vel into NE increment
    dvN = c*dvx - s*dvy
    dvE = s*dvx + c*dvy

    # discrete process model
    f = sp.Matrix([
        vN + dvN,
        vE + dvE,
        psi + dpsi,
        pN + dt*(vN + dvN),
        pE + dt*(vE + dvE),
    ])

    x = sp.Matrix([vN, vE, psi, pN, pE])
    w = sp.Matrix([dvx, dvy, dpsi])

    F = f.jacobian(x)
    G = f.jacobian(w)

    # Covariance prediction autocode
    # P is represented as symbolic symmetric matrix P00..P44
    P, _ = make_symm_P()
    dvxVar, dvyVar, dpsiVar = sp.symbols('dvxVar dvyVar dpsiVar', positive=True, real=True)
    Q = sp.diag(dvxVar, dvyVar, dpsiVar)
    P_pred = sp.simplify(F * P * F.T + G * Q * G.T)

    print('\n### P prediction (upper triangle):')
    pred_assign: list[tuple[str, sp.Expr]] = []
    for i in range(5):
        for j in range(i, 5):
            pred_assign.append((f'Pp[{i}][{j}]', P_pred[i, j]))
    # Full simplify was fine here (FLOW is where it tends to explode)
    emit_cse(pred_assign, temp_prefix='p', simplifier=sp.simplify)




    # =========================
    # Measurement update autocode
    # =========================

    # Position measurement update (2D):
    z_pos_n, z_pos_e = sp.symbols('z_pos_n z_pos_e', real=True)
    R_POS = sp.Symbol('R_POS', positive=True, real=True)  # (m)^2
    z_pos = sp.Matrix([z_pos_n, z_pos_e])


    # position measurement model
    h_pos = sp.Matrix([pN, pE])
    H_pos = h_pos.jacobian(x)
    
    """
    print('\n### H_pos (2x5):')
    for i in range(2):
        for j in range(5):
            print(f'H_pos[{i}][{j}] = {cexpr(H_pos[i,j])};')
    """
    innov_pos, S_pos, K_pos, x_pos, P_pos = ekf_update_2d(x, P, h_pos, z_pos, R_POS)

    print('\n### POS innov (pred - meas):')
    emit_cse([
        ('innov_pos[0]', innov_pos[0]),
        ('innov_pos[1]', innov_pos[1]),
    ], temp_prefix='pi', simplifier=sp.simplify)

    print('\n### POS S (2x2):')
    emit_cse([
        ('S_pos[0][0]', S_pos[0, 0]),
        ('S_pos[0][1]', S_pos[0, 1]),
        ('S_pos[1][0]', S_pos[1, 0]),
        ('S_pos[1][1]', S_pos[1, 1]),
    ], temp_prefix='ps', simplifier=sp.simplify)

    '''
    print('\n### POS K (5x2):')
    k_assign = []
    for i in range(5):
        for j in range(2):
            k_assign.append((f'K_pos[{i}][{j}]', K_pos[i, j]))
    emit_cse(k_assign, temp_prefix='pk')
    '''

    print('\n### POS x update:')
    emit_cse([
        ('xnew_pos[0]', x_pos[0]),
        ('xnew_pos[1]', x_pos[1]),
        ('xnew_pos[2]', x_pos[2]),
        ('xnew_pos[3]', x_pos[3]),
        ('xnew_pos[4]', x_pos[4]),
    ], temp_prefix='px', simplifier=sp.simplify)

    print('\n### POS P update (upper triangle):')
    p_assign = []
    for i in range(5):
        for j in range(i, 5):
            p_assign.append((f'Pnew_pos[{i}][{j}]', P_pos[i, j]))
    emit_cse(p_assign, temp_prefix='pp', simplifier=sp.simplify)


    # Flow measurement update (2D):
    # optical flow measurement model 
    # External parameters (not states)
    # we keep `rng` as an explicit input to keep the generated update compact,
    # but we also support computing it from main-EKF height and tilt:
    #   rng = hgt / (cos(roll)*cos(pitch))
    roll, pitch = sp.symbols('roll pitch', real=True)
    vD = sp.symbols('vD', real=True)
    p_rate, q_rate, r_rate = sp.symbols('p_rate q_rate r_rate', real=True)
    rx, ry, rz = sp.symbols('rx ry rz', real=True)
    rng = sp.symbols('rng', positive=True, real=True)

    # Height/tilt inputs and their variances (for augmenting S)
    hgt = sp.symbols('hgt', positive=True, real=True)  # height AGL (m)
    hgtVar = sp.Symbol('hgtVar', nonnegative=True, real=True)  # (m)^2
    rollVar = sp.Symbol('rollVar', nonnegative=True, real=True)  # (rad)^2
    pitchVar = sp.Symbol('pitchVar', nonnegative=True, real=True)  # (rad)^2

    # Body->NED rotation (321: yaw(Z), pitch(Y), roll(X))
    sr = sp.sin(roll)
    cr = sp.cos(roll)
    sp_ = sp.sin(pitch)
    cp = sp.cos(pitch)

    Rz = sp.Matrix([[c, -s, 0],
                    [s,  c, 0],
                    [0,  0, 1]])
    Ry = sp.Matrix([[cp, 0, sp_],
                    [0,  1, 0],
                    [-sp_, 0, cp]])
    Rx = sp.Matrix([[1, 0, 0],
                    [0, cr, -sr],
                    [0, sr,  cr]])

    # body->NED, then NED->body is transpose
    Rbn = Rz * Ry * Rx
    Tnb = Rbn.T

    velNED = sp.Matrix([vN, vE, vD])
    velBody = Tnb * velNED

    omega = sp.Matrix([p_rate, q_rate, r_rate])
    r_of = sp.Matrix([rx, ry, rz])
    omega_x_r = sp.Matrix([
        omega[1]*r_of[2] - omega[2]*r_of[1],
        omega[2]*r_of[0] - omega[0]*r_of[2],
        omega[0]*r_of[1] - omega[1]*r_of[0],
    ])

    relVelSensor = velBody + omega_x_r

    # EKF3 LOS prediction
    los_x = relVelSensor[1] / rng
    los_y = -relVelSensor[0] / rng
    h_flow = sp.Matrix([los_x, los_y])
    H_flow = h_flow.jacobian(x)

    # z_flow_x/y should be motion-compensated LOS rates: ofDataDelayed.flowRadXYcomp.(x,y)
    z_flow_x, z_flow_y = sp.symbols('z_flow_x z_flow_y', real=True)
    R_LOS = sp.Symbol('R_LOS', positive=True, real=True)  # (rad/s)^2

    z_flow = sp.Matrix([z_flow_x, z_flow_y])

    dlos_drng = sp.Matrix([sp.diff(h_flow[0], rng), sp.diff(h_flow[1], rng)])

    cosTilt = sp.cos(roll) * sp.cos(pitch)
    rng_from_hgt = hgt / cosTilt
    drng_dhgt = 1 / cosTilt
    drng_droll = sp.diff(rng_from_hgt, roll)
    drng_dpitch = sp.diff(rng_from_hgt, pitch)

    dlos_droll_partial = sp.Matrix([sp.diff(h_flow[0], roll), sp.diff(h_flow[1], roll)])
    dlos_dpitch_partial = sp.Matrix([sp.diff(h_flow[0], pitch), sp.diff(h_flow[1], pitch)])

    dlos_dhgt = dlos_drng * drng_dhgt
    dlos_droll = dlos_droll_partial + dlos_drng * drng_droll
    dlos_dpitch = dlos_dpitch_partial + dlos_drng * drng_dpitch

    S_add00_expr = dlos_dhgt[0]**2 * hgtVar + dlos_droll[0]**2 * rollVar + dlos_dpitch[0]**2 * pitchVar
    S_add01_expr = dlos_dhgt[0] * dlos_dhgt[1] * hgtVar + dlos_droll[0] * dlos_droll[1] * rollVar + dlos_dpitch[0] * dlos_dpitch[1] * pitchVar
    S_add11_expr = dlos_dhgt[1]**2 * hgtVar + dlos_droll[1]**2 * rollVar + dlos_dpitch[1]**2 * pitchVar

    # Keep the parameter-uncertainty terms symbolic in the EKF update to avoid huge expansion.
    S_add00, S_add01, S_add11 = sp.symbols('S_add00 S_add01 S_add11', real=True)
    extra_S_flow = sp.Matrix([[S_add00, S_add01], [S_add01, S_add11]])

    innov_flow, S_flow, K_flow, x_flow, P_flow = ekf_update_2d(x, P, h_flow, z_flow, R_LOS, extra_S=extra_S_flow)

    def flow_simplify(expr: sp.Expr) -> sp.Expr:
        # Keep flow expressions manageable by simplifying coefficients for the
        # linear terms in states/measurements, then do a cheap cleanup pass.
        expr = simplify_linear_coeffs(expr, vN, vE, z_flow_x, z_flow_y)
        return cheap_simplify(expr)

    print('\n### FLOW rng from height/tilt:')
    print('rng = hgt/(cos(roll)*cos(pitch));')

    print('\n### FLOW extra S from hgt/roll/pitch uncertainty:')
    emit_cse([
        ('S_add00', S_add00_expr),
        ('S_add01', S_add01_expr),
        ('S_add11', S_add11_expr),
    ], temp_prefix='fe', simplifier=flow_simplify)

    print('\n### FLOW innov (pred - meas):')
    emit_cse([
        ('innov_flow[0]', innov_flow[0]),
        ('innov_flow[1]', innov_flow[1]),
    ], temp_prefix='fi', simplifier=flow_simplify)

    print('\n### FLOW S (2x2):')
    emit_cse([
        ('S_flow[0][0]', S_flow[0, 0]),
        ('S_flow[0][1]', S_flow[0, 1]),
        ('S_flow[1][0]', S_flow[1, 0]),
        ('S_flow[1][1]', S_flow[1, 1]),
    ], temp_prefix='fs', simplifier=flow_simplify)

    print('\n### FLOW x update:')
    emit_cse([
        ('xnew_flow[0]', x_flow[0]),
        ('xnew_flow[1]', x_flow[1]),
        ('xnew_flow[2]', x_flow[2]),
        ('xnew_flow[3]', x_flow[3]),
        ('xnew_flow[4]', x_flow[4]),
    ], temp_prefix='fx', simplifier=flow_simplify)

    print('\n### FLOW P update (upper triangle):')
    p_assign: list[tuple[str, sp.Expr]] = []
    for i in range(5):
        for j in range(i, 5):
            p_assign.append((f'Pnew_flow[{i}][{j}]', P_flow[i, j]))
    emit_cse(p_assign, temp_prefix='fp', simplifier=flow_simplify)


if __name__ == '__main__':
    main()
