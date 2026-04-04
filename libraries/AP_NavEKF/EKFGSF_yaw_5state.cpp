/*
  Written from the ardupilot EKF3 implementation with 3 states

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <AP_HAL/AP_HAL.h>

#include "AP_NavEKF/EKFGSF_yaw_5state.h"
#include <AP_AHRS/AP_AHRS.h>

EKFGSF_yaw_5state::EKFGSF_yaw_5state()
{
    n_clips = 0;
}

ftype EKFGSF_yaw_5state::gaussianDensity(const uint8_t mdl_idx) const
{
    const ftype t2 = EKF[mdl_idx].S[0][0] * EKF[mdl_idx].S[1][1];
    const ftype t5 = EKF[mdl_idx].S[0][1] * EKF[mdl_idx].S[1][0];
    const ftype t3 = t2 - t5; // determinant
    const ftype t4 = 1.0f / MAX(t3, 1e-12f); // determinant inverse

    // inv(S)
    ftype invMat[2][2];
    invMat[0][0] =   t4 * EKF[mdl_idx].S[1][1];
    invMat[1][1] =   t4 * EKF[mdl_idx].S[0][0];
    invMat[0][1] = - t4 * EKF[mdl_idx].S[0][1];
    invMat[1][0] = - t4 * EKF[mdl_idx].S[1][0];

    // inv(S) * innovation
    ftype tempVec[2];
    tempVec[0] = invMat[0][0] * EKF[mdl_idx].innov[0] + invMat[0][1] * EKF[mdl_idx].innov[1];
    tempVec[1] = invMat[1][0] * EKF[mdl_idx].innov[0] + invMat[1][1] * EKF[mdl_idx].innov[1];

    // transpose(innovation) * inv(S) * innovation
    ftype normDist = tempVec[0] * EKF[mdl_idx].innov[0] + tempVec[1] * EKF[mdl_idx].innov[1];

    // convert from a normalised variance to a probability assuming a Gaussian distribution
    normDist = expf(-0.5f * normDist);
    normDist *= sqrtF(t4)/ M_2PI;
    return normDist;
}


void EKFGSF_yaw_5state::update(const Vector3F &delAng,
                        const Vector3F &delVel,
                        const ftype delAngDT,
                        const ftype delVelDT,
                        bool runEKF,
                        ftype TAS)
{
    // copy to class variables
    delta_angle = delAng;
    delta_velocity = delVel;
    angle_dt = delAngDT;
    velocity_dt = delVelDT;
    run_ekf_gsf = runEKF;
    true_airspeed = TAS;

    // Calculate a low pass filtered acceleration vector that will be used to keep the AHRS tilt aligned
    // The time constant of the filter is a fixed ratio relative to the time constant of the AHRS tilt correction loop
    const ftype filter_coef = fminF(EKFGSF_accelFiltRatio * delVelDT * EKFGSF_tiltGain, 1.0f);
    const Vector3F accel = delVel / fmaxF(delVelDT, 0.001f);
    ahrs_accel = ahrs_accel * (1.0f - filter_coef) + accel * filter_coef;

    // Iniitialise states and only when acceleration is close to 1g to prevent vehicle movement casuing a large initial tilt error
    if (!ahrs_tilt_aligned) {
        const ftype accel_norm_sq = accel.length_squared();
        const ftype upper_accel_limit = GRAVITY_MSS * 1.1f;
        const ftype lower_accel_limit = GRAVITY_MSS * 0.9f;
        const bool ok_to_align = ((accel_norm_sq > lower_accel_limit * lower_accel_limit &&
                                   accel_norm_sq < upper_accel_limit * upper_accel_limit));
        if (ok_to_align) {
            alignTilt();
            ahrs_tilt_aligned = true;
            ahrs_accel = accel;
        }
        return;
    }

    // Calculate common variables used by the AHRS prediction models
    ahrs_accel_norm = ahrs_accel.length();

    // Calculate AHRS acceleration fusion gain using a quadratic weighting function that is unity at 1g
    // and zero at the min and max g limits. This reduces the effect of large g transients on the attitude
    // esitmates.
    ftype EKFGSF_ahrs_ng = ahrs_accel_norm / GRAVITY_MSS;
    if (EKFGSF_ahrs_ng > 1.0f) {
        if (is_positive(true_airspeed)) {
            // When flying in fixed wing mode we need to allow for more positive g due to coordinated turns
            // Gain varies from unity at 1g to zero at 2g
            accel_gain = EKFGSF_tiltGain * sq(MAX(2.0f - EKFGSF_ahrs_ng, 0.0f));
        } else if (EKFGSF_ahrs_ng <= 1.5f) {
            // Gain varies from unity at 1g to zero at 1.5g
            accel_gain = EKFGSF_tiltGain * sq(3.0f - 2.0f * EKFGSF_ahrs_ng);
        } else {
            // Gain is zero above max g
            accel_gain = 0.0f;
        }
    } else if (EKFGSF_ahrs_ng > 0.5f) {
        // Gain varies from zero at 0.5g to unity at 1g
        accel_gain = EKFGSF_tiltGain * sq(2.0f * EKFGSF_ahrs_ng - 1.0f);
    } else {
        // Gain is zero below min g
        accel_gain = 0.0f;
    }

    // Always run the AHRS prediction cycle for each model
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        predict(mdl_idx);
    }

    if (vel_fuse_running && !run_ekf_gsf) {
        vel_fuse_running = false;
    }

    // Calculate a composite yaw as a weighted average of the states for each model.
    // To avoid issues with angle wrapping, the yaw state is converted to a vector with legnth
    // equal to the weighting value before it is summed.
    Vector2F yaw_vector;
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        yaw_vector[0] += GSF.weights[mdl_idx] * cosF(EKF[mdl_idx].X[2]);
        yaw_vector[1] += GSF.weights[mdl_idx] * sinF(EKF[mdl_idx].X[2]);
    }
    GSF.yaw = atan2F(yaw_vector[1],yaw_vector[0]);

    GSF.yaw_variance = 0.0f;
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        ftype yawDelta = wrap_PI(EKF[mdl_idx].X[2] - GSF.yaw);
        GSF.yaw_variance +=  GSF.weights[mdl_idx] * (EKF[mdl_idx].P[2][2] + sq(yawDelta));
    }
}


void EKFGSF_yaw_5state::fuseOFData(const Vector2F &flowRadXYcomp,
                            const ftype flowRadAcc,
                            const ftype heightAGL,
                            const ftype heightVar,
                            const ftype roll,
                            const ftype rollVar,
                            const ftype pitch,
                            const ftype pitchVar,
                            const ftype vD,
                            const Vector3F &bodyRates,
                            const Vector3F &sensorOffsetBody)
{
    // Compute slant range using main EKF height and tilt; do not depend on delayed rangefinder samples.
    // Assume flow sensor optical axis is aligned with body Z (down).
    const ftype cosTilt = cosF(roll) * cosF(pitch);
    if (!is_positive(heightAGL) || cosTilt < 0.1f) {
        return;
    }

    const ftype rng = heightAGL / cosTilt;
    // flowRadAcc is supplied by EKF3 parameter
    const ftype flowObsVar = sq(fmaxF(flowRadAcc, 1.0e-3f));

    if (run_ekf_gsf) {
        if (!vel_fuse_running) {
            resetEKFGSF();
            for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
                // Start from a neutral velocity hypothesis but with non-zero uncertainty.
                EKF[mdl_idx].X[0] = 0.0f;
                EKF[mdl_idx].X[1] = 0.0f;
                EKF[mdl_idx].P[0][0] = sq(5.0f);
                EKF[mdl_idx].P[1][1] = sq(5.0f);
            }
            alignYaw();
            vel_fuse_running = true;
        } else {
            ftype total_w = 0.0f;
            ftype newWeight[(uint8_t)N_MODELS_EKFGSF];
            bool state_update_failed = false;
            const ftype min_weight = 1e-5f;
            n_clips = 0;
            for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
                if (!correctOF(mdl_idx, flowRadXYcomp, flowObsVar, heightAGL, heightVar, roll, rollVar, pitch, pitchVar, vD, bodyRates, sensorOffsetBody, rng)) {
                    state_update_failed = true;
                }

                newWeight[mdl_idx] = gaussianDensity(mdl_idx) * GSF.weights[mdl_idx];
                if (newWeight[mdl_idx] < min_weight) {
                    n_clips++;
                    newWeight[mdl_idx] = min_weight;
                }
                total_w += newWeight[mdl_idx];
            }

            if (!state_update_failed) {
                if (vel_fuse_running && n_clips < N_MODELS_EKFGSF) {
                    ftype total_w_inv = 1.0f / total_w;
                    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
                        GSF.weights[mdl_idx] = newWeight[mdl_idx] * total_w_inv;
                    }
                } else {
                    resetEKFGSF();
                }
            }
        }
    }
}

void EKFGSF_yaw_5state::fusePosData(const Vector2F &pos, const ftype posAcc)
{
    const uint32_t now_ms = AP_HAL::millis();
    if (last_pos_fuse_ms > 0) {
        // Expected external-nav position rates are typically 8-30Hz.
        pos_meas_dt = fminF(fmaxF(1.0e-3f * (now_ms - last_pos_fuse_ms), 0.03f), 2.0f);
    } else {
        pos_meas_dt = 1.0f;
    }

    // enforce a velocity-equivalent floor like the 3-state GSF (0.5 m/s).
    const ftype pos_sigma_from_input = fmaxF(posAcc, 1.0e-3f);
    const ftype vel_sigma_floor = 0.5f;
    const ftype pos_sigma_floor = vel_sigma_floor * pos_meas_dt;
    const ftype posObsVar = sq(fmaxF(pos_sigma_from_input, pos_sigma_floor));

    if (run_ekf_gsf && vel_fuse_running) {
        bool state_update_failed = false;
        ftype total_w = 0.0f;
        ftype newWeight[(uint8_t)N_MODELS_EKFGSF];
        const ftype min_weight = 1e-5f;
        n_clips = 0;
        for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
            if (!correctPos(mdl_idx, pos, posObsVar)) {
                state_update_failed = true;
            }

            newWeight[mdl_idx] = gaussianDensity(mdl_idx) * GSF.weights[mdl_idx];
            if (newWeight[mdl_idx] < min_weight) {
                n_clips++;
                newWeight[mdl_idx] = min_weight;
            }
            total_w += newWeight[mdl_idx];
        }

        if (!state_update_failed) {
            if (vel_fuse_running && n_clips < N_MODELS_EKFGSF) {
                ftype total_w_inv = 1.0f / total_w;
                for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
                    GSF.weights[mdl_idx] = newWeight[mdl_idx] * total_w_inv;
                }
            } else {
                resetEKFGSF();
            }
            last_pos_fuse_ms = now_ms;
        }
    }
}

void EKFGSF_yaw_5state::predictAHRS(const uint8_t mdl_idx)
{
    // Generate attitude solution using simple complementary filter for the selected model

    // Calculate 'k' unit vector of earth frame rotated into body frame
    const Vector3F k{AHRS[mdl_idx].R[2][0], AHRS[mdl_idx].R[2][1], AHRS[mdl_idx].R[2][2]};

    // Calculate angular rate vector in rad/sec averaged across last sample interval
    const Vector3F ang_rate_delayed_raw { delta_angle / angle_dt };

    // Perform angular rate correction using accel data and reduce correction as accel magnitude moves away from 1 g (reduces drift when vehicle picked up and moved).
    // During fixed wing flight, compensate for centripetal acceleration assuming coordinated turns and X axis forward

    Vector3F tilt_error_gyro_correction; // (rad/sec)

    if (accel_gain > 0.0f) {

        Vector3F accel = ahrs_accel;

        if (is_positive(true_airspeed)) {
            // Calculate centripetal acceleration in body frame from cross product of body rate and body frame airspeed vector
            // NOTE: this assumes X axis is aligned with airspeed vector
            const Vector3F centripetal_accel_vec_bf {
                0.0f,
                ang_rate_delayed_raw[2] * true_airspeed,
                - ang_rate_delayed_raw[1] * true_airspeed
            };
            // Correct measured accel for centripetal acceleration
            accel -= centripetal_accel_vec_bf;
        }

        tilt_error_gyro_correction = (k % accel) * (accel_gain / ahrs_accel_norm);

    }

    // Gyro bias estimation
    const ftype gyro_bias_limit = radians(5.0f);
    const ftype spinRate_squared = ang_rate_delayed_raw.length_squared();
    if (spinRate_squared < sq(0.175f)) {
        AHRS[mdl_idx].gyro_bias -= tilt_error_gyro_correction * (EKFGSF_gyroBiasGain * angle_dt);

        // sanity check
        if (AHRS[mdl_idx].gyro_bias.is_nan()) {
            AHRS[mdl_idx].gyro_bias.zero();
        }

        for (uint8_t i = 0; i < 3; i++) {
            AHRS[mdl_idx].gyro_bias[i] = fminF(fmaxF(AHRS[mdl_idx].gyro_bias[i], -gyro_bias_limit), gyro_bias_limit);
        }
    }

    // Calculate the corrected body frame rotation vector for the last sample interval and apply to the rotation matrix
    const Vector3F ahrs_delta_angle = delta_angle + (tilt_error_gyro_correction - AHRS[mdl_idx].gyro_bias) * angle_dt;
    AHRS[mdl_idx].R = updateRotMat(AHRS[mdl_idx].R, ahrs_delta_angle);

}

void EKFGSF_yaw_5state::alignTilt()
{
    // Rotation matrix is constructed directly from acceleration measurement and will be the same for
    // all models so only need to calculate it once. Assumptions are:
    // 1) Yaw angle is zero - yaw is aligned later for each model when velocity fusion commences.
    // 2) The vehicle is not accelerating so all of the measured acceleration is due to gravity.

    // Calculate earth frame Down axis unit vector rotated into body frame
    Vector3F down_in_bf = -delta_velocity;
    down_in_bf.normalize();

    // Calculate earth frame North axis unit vector rotated into body frame, orthogonal to 'down_in_bf'
    // * operator is overloaded to provide a dot product
    const Vector3F i_vec_bf(1.0f,0.0f,0.0f);
    Vector3F north_in_bf = i_vec_bf - down_in_bf * (i_vec_bf * down_in_bf);
    north_in_bf.normalize();

    // Calculate earth frame East axis unit vector rotated into body frame, orthogonal to 'down_in_bf' and 'north_in_bf'
    // % operator is overloaded to provide a cross product
    const Vector3F east_in_bf = down_in_bf % north_in_bf;

    // Each column in a rotation matrix from earth frame to body frame represents the projection of the
    // corresponding earth frame unit vector rotated into the body frame, eg 'north_in_bf' would be the first column.
    // We need the rotation matrix from body frame to earth frame so the earth frame unit vectors rotated into body
    // frame are copied into corresponding rows instead to create the transpose.
    Matrix3F R;
    for (uint8_t col=0; col<3; col++) {
        R[0][col] = north_in_bf[col];
        R[1][col] = east_in_bf[col];
        R[2][col] = down_in_bf[col];
    }

    // record alignment
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        AHRS[mdl_idx].R = R;
        AHRS[mdl_idx].aligned = true;
    }
}

void EKFGSF_yaw_5state::alignYaw()
{
    // Align yaw angle for each model
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        if (fabsF(AHRS[mdl_idx].R[2][0]) < fabsF(AHRS[mdl_idx].R[2][1])) {
            // get the roll, pitch, yaw estimates from the rotation matrix using a  321 Tait-Bryan rotation sequence
            ftype roll,pitch,yaw;
            AHRS[mdl_idx].R.to_euler(&roll, &pitch, &yaw);

            // set the yaw angle
            yaw = wrap_PI(EKF[mdl_idx].X[2]);

            // update the body to earth frame rotation matrix
            AHRS[mdl_idx].R.from_euler(roll, pitch, yaw);

        } else {
            // Calculate the 312 Tait-Bryan rotation sequence that rotates from earth to body frame
            Vector3F euler312 = AHRS[mdl_idx].R.to_euler312();
            euler312[2] = wrap_PI(EKF[mdl_idx].X[2]); // first rotation (yaw) taken from EKF model state

            // update the body to earth frame rotation matrix
            AHRS[mdl_idx].R.from_euler312(euler312[0], euler312[1], euler312[2]);

        }
    }
}

// predict states and covariance for specified model index
void EKFGSF_yaw_5state::predict(const uint8_t mdl_idx)
{
    // generate an attitude reference using IMU data
    predictAHRS(mdl_idx);

    // we don't start running the EKF part of the algorithm until there are regular velocity observations
    if (!vel_fuse_running) {
        return;
    }

    // Calculate the yaw state using a projection onto the horizontal that avoids gimbal lock
    if (fabsF(AHRS[mdl_idx].R[2][0]) < fabsF(AHRS[mdl_idx].R[2][1])) {
        // use 321 Tait-Bryan rotation to define yaw state
        EKF[mdl_idx].X[2] = atan2F(AHRS[mdl_idx].R[1][0], AHRS[mdl_idx].R[0][0]);
    } else {
        // use 312 Tait-Bryan rotation to define yaw state
        EKF[mdl_idx].X[2] = atan2F(-AHRS[mdl_idx].R[0][1], AHRS[mdl_idx].R[1][1]); // first rotation (yaw)
    }

    // calculate delta velocity in a horizontal front-right frame
    const Vector3F del_vel_NED = AHRS[mdl_idx].R * delta_velocity;
    const ftype dvx =   del_vel_NED[0] * cosF(EKF[mdl_idx].X[2]) + del_vel_NED[1] * sinF(EKF[mdl_idx].X[2]);
    const ftype dvy = - del_vel_NED[0] * sinF(EKF[mdl_idx].X[2]) + del_vel_NED[1] * cosF(EKF[mdl_idx].X[2]);

    // sum delta velocities in earth frame:
    EKF[mdl_idx].X[0] += del_vel_NED[0];
    EKF[mdl_idx].X[1] += del_vel_NED[1];

    // update positions in earth frame via trapecoidal integration
    EKF[mdl_idx].X[3] += velocity_dt * EKF[mdl_idx].X[0];
    EKF[mdl_idx].X[4] += velocity_dt * EKF[mdl_idx].X[1];


    // predict covariance (5x5) - autocode from SymPy derivation

    // Local short variable name copies required for readability
    const ftype P00 = EKF[mdl_idx].P[0][0];
    const ftype P01 = EKF[mdl_idx].P[0][1];
    const ftype P02 = EKF[mdl_idx].P[0][2];
    const ftype P03 = EKF[mdl_idx].P[0][3];
    const ftype P04 = EKF[mdl_idx].P[0][4];
    const ftype P11 = EKF[mdl_idx].P[1][1];
    const ftype P12 = EKF[mdl_idx].P[1][2];
    const ftype P13 = EKF[mdl_idx].P[1][3];
    const ftype P14 = EKF[mdl_idx].P[1][4];
    const ftype P22 = EKF[mdl_idx].P[2][2];
    const ftype P23 = EKF[mdl_idx].P[2][3];
    const ftype P24 = EKF[mdl_idx].P[2][4];
    const ftype P33 = EKF[mdl_idx].P[3][3];
    const ftype P34 = EKF[mdl_idx].P[3][4];
    const ftype P44 = EKF[mdl_idx].P[4][4];

    const ftype dt = velocity_dt;

    // Use fixed values for delta velocity and delta angle process noise variances
    const ftype dvxVar = sq(EKFGSF_accelNoise * velocity_dt); // (m/s)^2
    const ftype dvyVar = dvxVar; // (m/s)^2
    const ftype dpsiVar = sq(EKFGSF_gyroNoise * angle_dt); // rad^2

    const ftype psi = EKF[mdl_idx].X[2];

    // ### P prediction (upper triangle):
    const ftype p0 = cosF(psi);
    const ftype p1 = sq(p0);
    const ftype p2 = dvxVar*p1;
    const ftype p3 = sinF(psi);
    const ftype p4 = sq(p3);
    const ftype p5 = dvyVar*p4;
    const ftype p6 = dvx*p3 + dvy*p0;
    const ftype p7 = P22*p6;
    const ftype p8 = P02 - p7;
    const ftype p9 = p6*p8;
    const ftype p10 = P02*p6;
    const ftype p11 = P00 - p10;
    const ftype p12 = dvx*p0 - dvy*p3;
    const ftype p13 = p12*p8;
    const ftype p14 = 0.5f*sinF(2.0f*psi);
    const ftype p15 = dvxVar*p14;
    const ftype p16 = dvyVar*p14;
    const ftype p17 = P12*p6;
    const ftype p18 = P01 - p17;
    const ftype p19 = P23*p6;
    const ftype p20 = P24*p6;
    const ftype p21 = dt*(p15 - p16);
    const ftype p22 = dvxVar*p4;
    const ftype p23 = dvyVar*p1;
    const ftype p24 = P22*p12;
    const ftype p25 = P12 + p24;
    const ftype p26 = p12*p25;
    const ftype p27 = P12*p12;
    const ftype p28 = P11 + p27;
    const ftype p29 = dt*p6;
    const ftype p30 = P24*p12;
    const ftype p31 = P02*dt + P23 - dt*p7;
    const ftype p32 = P12*dt + P24 + dt*p24;
    const ftype p33 = sq(dt);
    const ftype p34 = dt*p12;

    ftype Pp[5][5];
    Pp[0][0] = p11 + p2 + p5 - p9;
    Pp[0][1] = p13 + p15 - p16 + p18;
    Pp[0][2] = p8;
    Pp[0][3] = P03 + dt*p11 + dt*p2 + dt*p5 - dt*p9 - p19;
    Pp[0][4] = P04 + dt*p13 + dt*p18 - p20 + p21;
    Pp[1][1] = p22 + p23 + p26 + p28;
    Pp[1][2] = p25;
    Pp[1][3] = P13 + P23*p12 + dt*(P01 + P02*p12) + p21 - p25*p29;
    Pp[1][4] = P14 + dt*p22 + dt*p23 + dt*p26 + dt*p28 + p30;
    Pp[2][2] = P22 + dpsiVar;
    Pp[2][3] = p31;
    Pp[2][4] = p32;
    Pp[3][3] = P03*dt + P33 - dt*p19 + dt*(P00*dt + P03 - dt*p10) + p2*p33 - p29*p31 + p33*p5;
    Pp[3][4] = P04*dt + P34 - dt*p20 + dt*(P01*dt + P13 - dt*p17) + p15*p33 - p16*p33 + p31*p34;
    Pp[4][4] = P14*dt + P44 + dt*p30 + dt*(P11*dt + P14 + dt*p27) + p22*p33 + p23*p33 + p32*p34;

    const ftype min_var = 1e-6f;
    // Copy predicted covariance into P (fill both triangles)
    EKF[mdl_idx].P[0][0] = fmaxF(Pp[0][0], min_var);
    EKF[mdl_idx].P[0][1] = EKF[mdl_idx].P[1][0] = Pp[0][1];
    EKF[mdl_idx].P[0][2] = EKF[mdl_idx].P[2][0] = Pp[0][2];
    EKF[mdl_idx].P[0][3] = EKF[mdl_idx].P[3][0] = Pp[0][3];
    EKF[mdl_idx].P[0][4] = EKF[mdl_idx].P[4][0] = Pp[0][4];

    EKF[mdl_idx].P[1][1] = fmaxF(Pp[1][1], min_var);
    EKF[mdl_idx].P[1][2] = EKF[mdl_idx].P[2][1] = Pp[1][2];
    EKF[mdl_idx].P[1][3] = EKF[mdl_idx].P[3][1] = Pp[1][3];
    EKF[mdl_idx].P[1][4] = EKF[mdl_idx].P[4][1] = Pp[1][4];

    EKF[mdl_idx].P[2][2] = fmaxF(Pp[2][2], min_var);
    EKF[mdl_idx].P[2][3] = EKF[mdl_idx].P[3][2] = Pp[2][3];
    EKF[mdl_idx].P[2][4] = EKF[mdl_idx].P[4][2] = Pp[2][4];

    EKF[mdl_idx].P[3][3] = fmaxF(Pp[3][3], min_var);
    EKF[mdl_idx].P[3][4] = EKF[mdl_idx].P[4][3] = Pp[3][4];

    EKF[mdl_idx].P[4][4] = fmaxF(Pp[4][4], min_var);
}

bool EKFGSF_yaw_5state::correctOF(const uint8_t mdl_idx,
                           const Vector2F &flowRadXYcomp,
                           const ftype flowObsVar,
                           const ftype hgt,
                           const ftype hgtVar,
                           const ftype roll,
                           const ftype rollVar,
                           const ftype pitch,
                           const ftype pitchVar,
                           const ftype vD,
                           const Vector3F &bodyRates,
                           const Vector3F &sensorOffsetBody,
                           const ftype rng)
{
    // FLOW update autocode (SymPy-generated)
    if (!(rng > 0.1f)) {
        return false;
    }

    const ftype vN = EKF[mdl_idx].X[0];
    const ftype vE = EKF[mdl_idx].X[1];
    const ftype psi = EKF[mdl_idx].X[2];
    const ftype pN = EKF[mdl_idx].X[3];
    const ftype pE = EKF[mdl_idx].X[4];

    const ftype z_flow_x = flowRadXYcomp[0];
    const ftype z_flow_y = flowRadXYcomp[1];
    const ftype R_LOS = flowObsVar;

    const ftype p_rate = bodyRates[0];
    const ftype q_rate = bodyRates[1];
    const ftype r_rate = bodyRates[2];
    const ftype rx = sensorOffsetBody[0];
    const ftype ry = sensorOffsetBody[1];
    const ftype rz = sensorOffsetBody[2];

    // Local covariance copies
    const ftype P00 = EKF[mdl_idx].P[0][0];
    const ftype P01 = EKF[mdl_idx].P[0][1];
    const ftype P02 = EKF[mdl_idx].P[0][2];
    const ftype P03 = EKF[mdl_idx].P[0][3];
    const ftype P04 = EKF[mdl_idx].P[0][4];
    const ftype P11 = EKF[mdl_idx].P[1][1];
    const ftype P12 = EKF[mdl_idx].P[1][2];
    const ftype P13 = EKF[mdl_idx].P[1][3];
    const ftype P14 = EKF[mdl_idx].P[1][4];
    const ftype P22 = EKF[mdl_idx].P[2][2];
    const ftype P23 = EKF[mdl_idx].P[2][3];
    const ftype P24 = EKF[mdl_idx].P[2][4];
    const ftype P33 = EKF[mdl_idx].P[3][3];
    const ftype P34 = 0.5f * (EKF[mdl_idx].P[3][4] + EKF[mdl_idx].P[4][3]);
    const ftype P44 = EKF[mdl_idx].P[4][4];

    // ### FLOW rng from height/tilt:
    // rng = hgt/(cos(roll)*cos(pitch));

    // ### FLOW extra S from hgt/roll/pitch uncertainty:
    ftype S_add00;
    ftype S_add01;
    ftype S_add11;
    const ftype fe0 = sinF(roll);
    const ftype fe1 = cosF(pitch);
    const ftype fe2 = fe1*vD;
    const ftype fe3 = cosF(psi);
    const ftype fe4 = cosF(roll);
    const ftype fe5 = fe3*fe4;
    const ftype fe6 = sinF(pitch);
    const ftype fe7 = sinF(psi);
    const ftype fe8 = fe0*fe7;
    const ftype fe9 = fe4*fe7;
    const ftype fe10 = fe0*fe3;
    const ftype fe11 = fe0*fe2 - p_rate*rz + r_rate*rx + vE*(fe5 + fe6*fe8) + vN*(fe10*fe6 - fe9);
    const ftype fe12 = powf(fe4, -2.0f);
    const ftype fe13 = powf(fe1, -2.0f);
    const ftype fe14 = fe13/powf(rng, 4.0f);
    const ftype fe15 = fe12*fe14*hgtVar;
    const ftype fe16 = 1.0f/rng;
    const ftype fe17 = hgt/powf(rng, 2.0f);
    const ftype fe18 = fe11*fe17;
    const ftype fe19 = fe13*fe6/fe4;
    const ftype fe20 = -fe0*fe16*fe6*vD + fe1*fe10*fe16*vN + fe1*fe16*fe8*vE - fe18*fe19;
    const ftype fe21 = fe0*fe12/fe1;
    const ftype fe22 = fe16*fe2*fe4 + fe16*vE*(-fe10 + fe6*fe9) + fe16*vN*(fe5*fe6 + fe8) - fe18*fe21;
    const ftype fe23 = fe7*vE;
    const ftype fe24 = fe3*vN;
    const ftype fe25 = -fe1*fe23 - fe1*fe24 + fe6*vD - q_rate*rz + r_rate*ry;
    const ftype fe26 = fe17*fe25;
    const ftype fe27 = fe16*fe2 + fe16*fe23*fe6 + fe16*fe24*fe6 - fe19*fe26;
    const ftype fe28 = powf(fe25, 2.0f);
    S_add00 = powf(fe11, 2.0f)*fe15 + powf(fe20, 2.0f)*pitchVar + powf(fe22, 2.0f)*rollVar;
    S_add01 = fe11*fe15*fe25 + fe20*fe27*pitchVar - fe21*fe22*fe26*rollVar;
    S_add11 = (powf(fe0, 2.0f)*fe14*fe28*powf(hgt, 2.0f)*rollVar + fe15*fe28*powf(fe4, 4.0f) + powf(fe27, 2.0f)*powf(fe4, 4.0f)*pitchVar)/powf(fe4, 4.0f);

    // ### FLOW innov (pred - meas):
    const ftype fi0 = 1.0f/rng;
    const ftype fi1 = sinF(roll);
    const ftype fi2 = cosF(pitch);
    const ftype fi3 = cosF(psi);
    const ftype fi4 = cosF(roll);
    const ftype fi5 = sinF(psi);
    const ftype fi6 = sinF(pitch);
    const ftype fi7 = fi1*fi6;
    EKF[mdl_idx].innov_flow[0] = fi0*vE*(fi3*fi4 + fi5*fi7) + fi0*vN*(fi3*fi7 - fi4*fi5) + fi0*(fi1*fi2*vD - p_rate*rz + r_rate*rx) - z_flow_x;
    EKF[mdl_idx].innov_flow[1] = -fi0*fi2*fi3*vN - fi0*fi2*fi5*vE + fi0*(fi6*vD - q_rate*rz + r_rate*ry) - z_flow_y;
    // Convert LOS-rate innovation to a horizontal velocity proxy (m/s) for quality gating.
    EKF[mdl_idx].innov_vel[0] = EKF[mdl_idx].innov_flow[0] * rng;
    EKF[mdl_idx].innov_vel[1] = EKF[mdl_idx].innov_flow[1] * rng;
    EKF[mdl_idx].innov_flow_vel[0] = EKF[mdl_idx].innov_vel[0];
    EKF[mdl_idx].innov_flow_vel[1] = EKF[mdl_idx].innov_vel[1];
    EKF[mdl_idx].innov[0] = EKF[mdl_idx].innov_flow[0];
    EKF[mdl_idx].innov[1] = EKF[mdl_idx].innov_flow[1];

    // ### FLOW S (2x2):
    const ftype fs0 = 1.0f/rng;
    const ftype fs1 = cosF(psi);
    const ftype fs2 = cosF(roll);
    const ftype fs3 = sinF(psi);
    const ftype fs4 = sinF(pitch)*sinF(roll);
    const ftype fs5 = fs1*fs2 + fs3*fs4;
    const ftype fs6 = fs0*fs5;
    const ftype fs7 = fs1*fs4 - fs2*fs3;
    const ftype fs8 = fs0*fs7;
    const ftype fs9 = fs0*(-fs5*vN + fs7*vE);
    const ftype fs10 = P01*fs8 + P11*fs6 + P12*fs9;
    const ftype fs11 = P00*fs8 + P01*fs6 + P02*fs9;
    const ftype fs12 = P02*fs8 + P12*fs6 + P22*fs9;
    const ftype fs13 = cosF(pitch);
    const ftype fs14 = fs1*fs13;
    const ftype fs15 = fs0*fs14;
    const ftype fs16 = fs0*fs13*fs3;
    const ftype fs17 = fs13*fs3*vN - fs14*vE;
    const ftype fs18 = -P01*fs15 - P11*fs16 + P12*fs0*fs17;
    const ftype fs19 = -P00*fs15 - P01*fs16 + P02*fs0*fs17;
    const ftype fs20 = -P02*fs15 - P12*fs16 + P22*fs0*fs17;
    const ftype fs21 = fs0*fs17;
    EKF[mdl_idx].S_flow[0][0] = R_LOS + S_add00 + fs10*fs6 + fs11*fs8 + fs12*fs9;
    EKF[mdl_idx].S_flow[0][1] = S_add01 + fs18*fs6 + fs19*fs8 + fs20*fs9;
    EKF[mdl_idx].S_flow[1][0] = S_add01 - fs10*fs16 - fs11*fs15 + fs12*fs21;
    EKF[mdl_idx].S_flow[1][1] = R_LOS + S_add11 - fs15*fs19 - fs16*fs18 + fs20*fs21;
    EKF[mdl_idx].S[0][0] = EKF[mdl_idx].S_flow[0][0];
    EKF[mdl_idx].S[0][1] = EKF[mdl_idx].S_flow[0][1];
    EKF[mdl_idx].S[1][0] = EKF[mdl_idx].S_flow[1][0];
    EKF[mdl_idx].S[1][1] = EKF[mdl_idx].S_flow[1][1];

    // chi-square innovation compression: clip effective correction to 5-sigma.
    ftype innov_comp_scale_factor = 1.0f;
    ftype S_det_inv = EKF[mdl_idx].S[0][0]*EKF[mdl_idx].S[1][1] - EKF[mdl_idx].S[0][1]*EKF[mdl_idx].S[1][0];
    if (fabsF(S_det_inv) > 1E-6f) {
        S_det_inv = 1.0f / S_det_inv;
        const ftype S_inv_NN = EKF[mdl_idx].S[1][1] * S_det_inv;
        const ftype S_inv_EE = EKF[mdl_idx].S[0][0] * S_det_inv;
        const ftype S_inv_NE = EKF[mdl_idx].S[0][1] * S_det_inv;
        const ftype test_ratio = EKF[mdl_idx].innov[0]*(EKF[mdl_idx].innov[0]*S_inv_NN + EKF[mdl_idx].innov[1]*S_inv_NE) +
                                 EKF[mdl_idx].innov[1]*(EKF[mdl_idx].innov[0]*S_inv_NE + EKF[mdl_idx].innov[1]*S_inv_EE);
        if (test_ratio > 25.0f) {
            innov_comp_scale_factor = sqrtF(25.0f / test_ratio);
        }
    }
    const ftype z_flow_x_compressed = z_flow_x + (1.0f - innov_comp_scale_factor) * EKF[mdl_idx].innov[0];
    const ftype z_flow_y_compressed = z_flow_y + (1.0f - innov_comp_scale_factor) * EKF[mdl_idx].innov[1];

    // ### FLOW x update:
    const ftype fx0 = 1.0f/rng;
    const ftype fx1 = sinF(roll);
    const ftype fx2 = cosF(pitch);
    const ftype fx3 = cosF(psi);
    const ftype fx4 = cosF(roll);
    const ftype fx5 = sinF(psi);
    const ftype fx6 = sinF(pitch);
    const ftype fx7 = fx1*fx6;
    const ftype fx8 = fx3*fx4 + fx5*fx7;
    const ftype fx9 = fx3*fx7 - fx4*fx5;
    const ftype fx10 = fx0*fx8*vE + fx0*fx9*vN + fx0*(fx1*fx2*vD - p_rate*rz + r_rate*rx) - z_flow_x_compressed;
    const ftype fx11 = fx2*fx3;
    const ftype fx12 = fx0*fx11;
    const ftype fx13 = fx2*fx5;
    const ftype fx14 = fx0*fx13;
    const ftype fx15 = -fx11*vE + fx2*fx5*vN;
    const ftype fx16 = -P00*fx12 - P01*fx14 + P02*fx0*fx15;
    const ftype fx17 = -P01*fx12 - P11*fx14 + P12*fx0*fx15;
    const ftype fx18 = -P02*fx12 - P12*fx14 + P22*fx0*fx15;
    const ftype fx19 = fx0*fx15;
    const ftype fx20 = R_LOS + S_add11 - fx12*fx16 - fx14*fx17 + fx18*fx19;
    const ftype fx21 = fx0*fx8;
    const ftype fx22 = fx0*fx9;
    const ftype fx23 = fx0*(-fx8*vN + fx9*vE);
    const ftype fx24 = P00*fx22 + P01*fx21 + P02*fx23;
    const ftype fx25 = S_add01 + fx16*fx22 + fx17*fx21 + fx18*fx23;
    const ftype fx26 = P01*fx22 + P11*fx21 + P12*fx23;
    const ftype fx27 = P02*fx22 + P12*fx21 + P22*fx23;
    const ftype fx28 = S_add01 - fx12*fx24 - fx14*fx26 + fx19*fx27;
    const ftype fx29 = R_LOS + S_add00 + fx21*fx26 + fx22*fx24 + fx23*fx27;
    const ftype fx30_den = fx20*fx29 - fx25*fx28;
    if (fabsF(fx30_den) < 1e-12f) {
        return false;
    }
    const ftype fx30 = 1.0f/fx30_den;
    const ftype fx31 = fx24*fx30;
    const ftype fx32 = -fx28;
    const ftype fx33 = fx16*fx30;
    const ftype fx34 = -fx0*fx11*vN - fx0*fx13*vE + fx0*(fx6*vD - q_rate*rz + r_rate*ry) - z_flow_y_compressed;
    const ftype fx35 = -fx25;
    const ftype fx36 = fx26*fx30;
    const ftype fx37 = fx17*fx30;
    const ftype fx38 = fx27*fx30;
    const ftype fx39 = fx18*fx30;
    const ftype fx40 = fx30*(P03*fx22 + P13*fx21 + P23*fx23);
    const ftype fx41 = fx30*(-P03*fx12 - P13*fx14 + P23*fx0*fx15);
    const ftype fx42 = fx30*(P04*fx22 + P14*fx21 + P24*fx23);
    const ftype fx43 = fx30*(-P04*fx12 - P14*fx14 + P24*fx0*fx15);
    const ftype xnew0 = -fx10*fx20*fx31 - fx10*fx32*fx33 - fx29*fx33*fx34 - fx31*fx34*fx35 + vN;
    const ftype xnew1 = -fx10*fx20*fx36 - fx10*fx32*fx37 - fx29*fx34*fx37 - fx34*fx35*fx36 + vE;
    const ftype xnew2 = -fx10*fx20*fx38 - fx10*fx32*fx39 - fx29*fx34*fx39 - fx34*fx35*fx38 + psi;
    const ftype xnew3 = -fx10*fx20*fx40 - fx10*fx32*fx41 - fx29*fx34*fx41 - fx34*fx35*fx40 + pN;
    const ftype xnew4 = -fx10*fx20*fx42 - fx10*fx32*fx43 - fx29*fx34*fx43 - fx34*fx35*fx42 + pE;

    const ftype yaw_prev = psi;
    const ftype yaw_new = wrap_PI(xnew2);
    const ftype yaw_delta = wrap_PI(yaw_new - yaw_prev);

    EKF[mdl_idx].X[0] = xnew0;
    EKF[mdl_idx].X[1] = xnew1;
    EKF[mdl_idx].X[2] = yaw_new;
    EKF[mdl_idx].X[3] = xnew3;
    EKF[mdl_idx].X[4] = xnew4;

    // apply the change in yaw angle to the AHRS taking advantage of sparseness in the yaw rotation matrix
    const ftype cos_yaw = cosF(yaw_delta);
    const ftype sin_yaw = sinF(yaw_delta);
    ftype R_prev[2][3];
    memcpy(&R_prev, &AHRS[mdl_idx].R, sizeof(R_prev));
    AHRS[mdl_idx].R[0][0] = R_prev[0][0] * cos_yaw - R_prev[1][0] * sin_yaw;
    AHRS[mdl_idx].R[0][1] = R_prev[0][1] * cos_yaw - R_prev[1][1] * sin_yaw;
    AHRS[mdl_idx].R[0][2] = R_prev[0][2] * cos_yaw - R_prev[1][2] * sin_yaw;
    AHRS[mdl_idx].R[1][0] = R_prev[0][0] * sin_yaw + R_prev[1][0] * cos_yaw;
    AHRS[mdl_idx].R[1][1] = R_prev[0][1] * sin_yaw + R_prev[1][1] * cos_yaw;
    AHRS[mdl_idx].R[1][2] = R_prev[0][2] * sin_yaw + R_prev[1][2] * cos_yaw;

    // ### FLOW P update (upper triangle):
    const ftype fp0 = 1.0f/rng;
    const ftype fp1 = cosF(psi);
    const ftype fp2 = cosF(roll);
    const ftype fp3 = sinF(psi);
    const ftype fp4 = sinF(pitch)*sinF(roll);
    const ftype fp5 = fp1*fp2 + fp3*fp4;
    const ftype fp6 = fp0*fp5;
    const ftype fp7 = fp1*fp4 - fp2*fp3;
    const ftype fp8 = fp0*fp7;
    const ftype fp9 = fp0*(-fp5*vN + fp7*vE);
    const ftype fp10 = P00*fp8 + P01*fp6 + P02*fp9;
    const ftype fp11 = cosF(pitch);
    const ftype fp12 = fp1*fp11;
    const ftype fp13 = fp0*fp12;
    const ftype fp14 = fp0*fp11*fp3;
    const ftype fp15 = fp11*fp3*vN - fp12*vE;
    const ftype fp16 = -P00*fp13 - P01*fp14 + P02*fp0*fp15;
    const ftype fp17 = -P01*fp13 - P11*fp14 + P12*fp0*fp15;
    const ftype fp18 = -P02*fp13 - P12*fp14 + P22*fp0*fp15;
    const ftype fp19 = fp0*fp15;
    const ftype fp20 = R_LOS + S_add11 - fp13*fp16 - fp14*fp17 + fp18*fp19;
    const ftype fp21 = S_add01 + fp16*fp8 + fp17*fp6 + fp18*fp9;
    const ftype fp22 = P01*fp8 + P11*fp6 + P12*fp9;
    const ftype fp23 = P02*fp8 + P12*fp6 + P22*fp9;
    const ftype fp24 = S_add01 - fp10*fp13 - fp14*fp22 + fp19*fp23;
    const ftype fp25 = R_LOS + S_add00 + fp10*fp8 + fp22*fp6 + fp23*fp9;
    const ftype fp26_den = fp20*fp25 - fp21*fp24;
    if (fabsF(fp26_den) < 1e-12f) {
        return false;
    }
    const ftype fp26 = 1.0f/fp26_den;
    const ftype fp27 = fp10*fp26;
    const ftype fp28 = -fp24;
    const ftype fp29 = fp16*fp26;
    const ftype fp30 = fp20*fp27 + fp28*fp29;
    const ftype fp31 = -fp21;
    const ftype fp32 = fp25*fp29 + fp27*fp31;
    const ftype fp33 = P03*fp8 + P13*fp6 + P23*fp9;
    const ftype fp34 = -P03*fp13 - P13*fp14 + P23*fp0*fp15;
    const ftype fp35 = P04*fp8 + P14*fp6 + P24*fp9;
    const ftype fp36 = -P04*fp13 - P14*fp14 + P24*fp0*fp15;
    const ftype fp37 = fp22*fp26;
    const ftype fp38 = fp17*fp26;
    const ftype fp39 = fp20*fp37 + fp28*fp38;
    const ftype fp40 = fp25*fp38 + fp31*fp37;
    const ftype fp41 = fp23*fp26;
    const ftype fp42 = fp18*fp26;
    const ftype fp43 = fp20*fp41 + fp28*fp42;
    const ftype fp44 = fp25*fp42 + fp31*fp41;
    const ftype fp45 = fp26*fp33;
    const ftype fp46 = fp26*fp34;
    const ftype fp47 = fp20*fp45 + fp28*fp46;
    const ftype fp48 = fp25*fp46 + fp31*fp45;
    const ftype fp49 = fp26*fp35;
    const ftype fp50 = fp26*fp36;

    ftype Pnew_flow[5][5];
    Pnew_flow[0][0] = P00 - fp10*fp30 - fp16*fp32;
    Pnew_flow[0][1] = P01 - fp17*fp32 - fp22*fp30;
    Pnew_flow[0][2] = P02 - fp18*fp32 - fp23*fp30;
    Pnew_flow[0][3] = P03 - fp30*fp33 - fp32*fp34;
    Pnew_flow[0][4] = P04 - fp30*fp35 - fp32*fp36;
    Pnew_flow[1][1] = P11 - fp17*fp40 - fp22*fp39;
    Pnew_flow[1][2] = P12 - fp18*fp40 - fp23*fp39;
    Pnew_flow[1][3] = P13 - fp33*fp39 - fp34*fp40;
    Pnew_flow[1][4] = P14 - fp35*fp39 - fp36*fp40;
    Pnew_flow[2][2] = P22 - fp18*fp44 - fp23*fp43;
    Pnew_flow[2][3] = P23 - fp33*fp43 - fp34*fp44;
    Pnew_flow[2][4] = P24 - fp35*fp43 - fp36*fp44;
    Pnew_flow[3][3] = P33 - fp33*fp47 - fp34*fp48;
    Pnew_flow[3][4] = P34 - fp35*fp47 - fp36*fp48;
    Pnew_flow[4][4] = P44 - fp20*fp35*fp49 - fp25*fp36*fp50 - fp28*fp35*fp50 - fp31*fp36*fp49;

    // Write back covariance (fill symmetric elements)
    EKF[mdl_idx].P[0][0] = Pnew_flow[0][0];
    EKF[mdl_idx].P[0][1] = EKF[mdl_idx].P[1][0] = Pnew_flow[0][1];
    EKF[mdl_idx].P[0][2] = EKF[mdl_idx].P[2][0] = Pnew_flow[0][2];
    EKF[mdl_idx].P[0][3] = EKF[mdl_idx].P[3][0] = Pnew_flow[0][3];
    EKF[mdl_idx].P[0][4] = EKF[mdl_idx].P[4][0] = Pnew_flow[0][4];

    EKF[mdl_idx].P[1][1] = Pnew_flow[1][1];
    EKF[mdl_idx].P[1][2] = EKF[mdl_idx].P[2][1] = Pnew_flow[1][2];
    EKF[mdl_idx].P[1][3] = EKF[mdl_idx].P[3][1] = Pnew_flow[1][3];
    EKF[mdl_idx].P[1][4] = EKF[mdl_idx].P[4][1] = Pnew_flow[1][4];

    EKF[mdl_idx].P[2][2] = Pnew_flow[2][2];
    EKF[mdl_idx].P[2][3] = EKF[mdl_idx].P[3][2] = Pnew_flow[2][3];
    EKF[mdl_idx].P[2][4] = EKF[mdl_idx].P[4][2] = Pnew_flow[2][4];

    EKF[mdl_idx].P[3][3] = Pnew_flow[3][3];
    EKF[mdl_idx].P[3][4] = EKF[mdl_idx].P[4][3] = Pnew_flow[3][4];

    EKF[mdl_idx].P[4][4] = Pnew_flow[4][4];

    const ftype min_var = 1e-6f;
    for (uint8_t i = 0; i < 5; i++) {
        EKF[mdl_idx].P[i][i] = fmaxF(EKF[mdl_idx].P[i][i], min_var);
    }

    return true;
}

bool EKFGSF_yaw_5state::correctPos(const uint8_t mdl_idx, const Vector2F &pos, const ftype posObsVar)
{
    // POS update autocode (SymPy-generated)
    const ftype vN = EKF[mdl_idx].X[0];
    const ftype vE = EKF[mdl_idx].X[1];
    const ftype psi = EKF[mdl_idx].X[2];
    const ftype pN = EKF[mdl_idx].X[3];
    const ftype pE = EKF[mdl_idx].X[4];

    const ftype z_pos_n = pos[0];
    const ftype z_pos_e = pos[1];
    const ftype R_POS = posObsVar;

    // Local covariance copies
    const ftype P00 = EKF[mdl_idx].P[0][0];
    const ftype P01 = EKF[mdl_idx].P[0][1];
    const ftype P02 = EKF[mdl_idx].P[0][2];
    const ftype P03 = EKF[mdl_idx].P[0][3];
    const ftype P04 = EKF[mdl_idx].P[0][4];
    const ftype P11 = EKF[mdl_idx].P[1][1];
    const ftype P12 = EKF[mdl_idx].P[1][2];
    const ftype P13 = EKF[mdl_idx].P[1][3];
    const ftype P14 = EKF[mdl_idx].P[1][4];
    const ftype P22 = EKF[mdl_idx].P[2][2];
    const ftype P23 = EKF[mdl_idx].P[2][3];
    const ftype P24 = EKF[mdl_idx].P[2][4];
    const ftype P33 = EKF[mdl_idx].P[3][3];
    const ftype P34 = 0.5f * (EKF[mdl_idx].P[3][4] + EKF[mdl_idx].P[4][3]);
    const ftype P44 = EKF[mdl_idx].P[4][4];

    // ### POS innov (pred - meas):
    EKF[mdl_idx].innov_pos[0] = pN - z_pos_n;
    EKF[mdl_idx].innov_pos[1] = pE - z_pos_e;
    // Convert position innovation to a velocity proxy (m/s) using the measured position sample interval.
    EKF[mdl_idx].innov_vel[0] = EKF[mdl_idx].innov_pos[0] / pos_meas_dt;
    EKF[mdl_idx].innov_vel[1] = EKF[mdl_idx].innov_pos[1] / pos_meas_dt;
    EKF[mdl_idx].innov[0] = EKF[mdl_idx].innov_pos[0];
    EKF[mdl_idx].innov[1] = EKF[mdl_idx].innov_pos[1];

    // ### POS S (2x2):
    EKF[mdl_idx].S_pos[0][0] = P33 + R_POS;
    EKF[mdl_idx].S_pos[0][1] = P34;
    EKF[mdl_idx].S_pos[1][0] = P34;
    EKF[mdl_idx].S_pos[1][1] = P44 + R_POS;
    EKF[mdl_idx].S[0][0] = EKF[mdl_idx].S_pos[0][0];
    EKF[mdl_idx].S[0][1] = EKF[mdl_idx].S_pos[0][1];
    EKF[mdl_idx].S[1][0] = EKF[mdl_idx].S_pos[1][0];
    EKF[mdl_idx].S[1][1] = EKF[mdl_idx].S_pos[1][1];

    // Legacy-style chi-square innovation compression: clip effective correction to 5-sigma.
    // This scales state correction only; covariance update remains unchanged.
    ftype innov_comp_scale_factor = 1.0f;
    ftype S_det_inv = EKF[mdl_idx].S[0][0]*EKF[mdl_idx].S[1][1] - EKF[mdl_idx].S[0][1]*EKF[mdl_idx].S[1][0];
    if (fabsF(S_det_inv) > 1E-6f) {
        S_det_inv = 1.0f / S_det_inv;
        const ftype S_inv_NN = EKF[mdl_idx].S[1][1] * S_det_inv;
        const ftype S_inv_EE = EKF[mdl_idx].S[0][0] * S_det_inv;
        const ftype S_inv_NE = EKF[mdl_idx].S[0][1] * S_det_inv;
        const ftype test_ratio = EKF[mdl_idx].innov[0]*(EKF[mdl_idx].innov[0]*S_inv_NN + EKF[mdl_idx].innov[1]*S_inv_NE) +
                                 EKF[mdl_idx].innov[1]*(EKF[mdl_idx].innov[0]*S_inv_NE + EKF[mdl_idx].innov[1]*S_inv_EE);
        if (test_ratio > 25.0f) {
            innov_comp_scale_factor = sqrtF(25.0f / test_ratio);
        }
    }
    const ftype z_pos_n_compressed = z_pos_n + (1.0f - innov_comp_scale_factor) * EKF[mdl_idx].innov[0];
    const ftype z_pos_e_compressed = z_pos_e + (1.0f - innov_comp_scale_factor) * EKF[mdl_idx].innov[1];

    // Guard against bad conditioning (det(S) too small)
    const ftype detS = (P33 + R_POS) * (P44 + R_POS) - sq(P34);
    if (fabsF(detS) < 1e-12f) {
        return false;
    }

    // ### POS x update:
    const ftype px0 = pE - z_pos_e_compressed;
    const ftype px1 = sq(P34);
    const ftype px2 = 1.0f/(P33*P44 + P33*R_POS + P44*R_POS + sq(R_POS) - px1);
    const ftype px3 = P03*px2;
    const ftype px4 = P33 + R_POS;
    const ftype px5 = pN - z_pos_n_compressed;
    const ftype px6 = P44 + R_POS;
    const ftype px7 = P13*px2;
    const ftype px8 = P23*px2;
    const ftype px9 = P34*px2;
    const ftype px10 = px1*px2;
    const ftype xnew0 = -px0*(P04*px2*px4 - P34*px3) + px5*(P04*P34*px2 - px3*px6) + vN;
    const ftype xnew1 = -px0*(P14*px2*px4 - P34*px7) + px5*(P14*P34*px2 - px6*px7) + vE;
    const ftype xnew2 = psi - px0*(P24*px2*px4 - P34*px8) + px5*(P24*P34*px2 - px6*px8);
    const ftype xnew3 = pN + px0*(P33*px9 - P34*px2*px4) - px5*(P33*px2*px6 - px10);
    const ftype xnew4 = pE - px0*(P44*px2*px4 - px10) - px5*(P34*px2*px6 - P44*px9);

    const ftype yaw_prev = psi;
    const ftype yaw_new = wrap_PI(xnew2);
    const ftype yaw_delta = wrap_PI(yaw_new - yaw_prev);

    EKF[mdl_idx].X[0] = xnew0;
    EKF[mdl_idx].X[1] = xnew1;
    EKF[mdl_idx].X[2] = yaw_new;
    EKF[mdl_idx].X[3] = xnew3;
    EKF[mdl_idx].X[4] = xnew4;

    // apply the change in yaw angle to the AHRS taking advantage of sparseness in the yaw rotation matrix
    const ftype cos_yaw = cosF(yaw_delta);
    const ftype sin_yaw = sinF(yaw_delta);
    ftype R_prev[2][3];
    memcpy(&R_prev, &AHRS[mdl_idx].R, sizeof(R_prev));
    AHRS[mdl_idx].R[0][0] = R_prev[0][0] * cos_yaw - R_prev[1][0] * sin_yaw;
    AHRS[mdl_idx].R[0][1] = R_prev[0][1] * cos_yaw - R_prev[1][1] * sin_yaw;
    AHRS[mdl_idx].R[0][2] = R_prev[0][2] * cos_yaw - R_prev[1][2] * sin_yaw;
    AHRS[mdl_idx].R[1][0] = R_prev[0][0] * sin_yaw + R_prev[1][0] * cos_yaw;
    AHRS[mdl_idx].R[1][1] = R_prev[0][1] * sin_yaw + R_prev[1][1] * cos_yaw;
    AHRS[mdl_idx].R[1][2] = R_prev[0][2] * sin_yaw + R_prev[1][2] * cos_yaw;

    // ### POS P update (upper triangle):
    const ftype pp0 = sq(P34);
    const ftype pp1 = 1.0f/(P33*P44 + P33*R_POS + P44*R_POS + sq(R_POS) - pp0);
    const ftype pp2 = P44 + R_POS;
    const ftype pp3 = P03*pp1;
    const ftype pp4 = -P04*P34*pp1 + pp2*pp3;
    const ftype pp5 = P33 + R_POS;
    const ftype pp6 = P04*pp1*pp5 - P34*pp3;
    const ftype pp7 = P13*pp1;
    const ftype pp8 = -P14*P34*pp1 + pp2*pp7;
    const ftype pp9 = P14*pp1*pp5 - P34*pp7;
    const ftype pp10 = P23*pp1;
    const ftype pp11 = -P24*P34*pp1 + pp10*pp2;
    const ftype pp12 = P24*pp1*pp5 - P34*pp10;
    const ftype pp13 = P33*pp1;
    const ftype pp14 = P34*(pp1*pp5 - pp13);
    const ftype pp15 = pp0*pp1;
    const ftype pp16 = pp13*pp2 - pp15;

    ftype Pnew_pos[5][5];
    Pnew_pos[0][0] = P00 - P03*pp4 - P04*pp6;
    Pnew_pos[0][1] = P01 - P13*pp4 - P14*pp6;
    Pnew_pos[0][2] = P02 - P23*pp4 - P24*pp6;
    Pnew_pos[0][3] = P03 - P33*pp4 - P34*pp6;
    Pnew_pos[0][4] = P04 - P34*pp4 - P44*pp6;
    Pnew_pos[1][1] = P11 - P13*pp8 - P14*pp9;
    Pnew_pos[1][2] = P12 - P23*pp8 - P24*pp9;
    Pnew_pos[1][3] = P13 - P33*pp8 - P34*pp9;
    Pnew_pos[1][4] = P14 - P34*pp8 - P44*pp9;
    Pnew_pos[2][2] = P22 - P23*pp11 - P24*pp12;
    Pnew_pos[2][3] = P23 - P33*pp11 - P34*pp12;
    Pnew_pos[2][4] = P24 - P34*pp11 - P44*pp12;
    Pnew_pos[3][3] = -P33*pp16 + P33 - P34*pp14;
    Pnew_pos[3][4] = -P34*pp16 + P34 - P44*pp14;
    Pnew_pos[4][4] = sq(P34)*pp1*(P44 - pp2) - P44*(P44*pp1*pp5 - pp15) + P44;

    // Write back covariance (fill symmetric elements)
    EKF[mdl_idx].P[0][0] = Pnew_pos[0][0];
    EKF[mdl_idx].P[0][1] = EKF[mdl_idx].P[1][0] = Pnew_pos[0][1];
    EKF[mdl_idx].P[0][2] = EKF[mdl_idx].P[2][0] = Pnew_pos[0][2];
    EKF[mdl_idx].P[0][3] = EKF[mdl_idx].P[3][0] = Pnew_pos[0][3];
    EKF[mdl_idx].P[0][4] = EKF[mdl_idx].P[4][0] = Pnew_pos[0][4];
    EKF[mdl_idx].P[1][1] = Pnew_pos[1][1];
    EKF[mdl_idx].P[1][2] = EKF[mdl_idx].P[2][1] = Pnew_pos[1][2];
    EKF[mdl_idx].P[1][3] = EKF[mdl_idx].P[3][1] = Pnew_pos[1][3];
    EKF[mdl_idx].P[1][4] = EKF[mdl_idx].P[4][1] = Pnew_pos[1][4];
    EKF[mdl_idx].P[2][2] = Pnew_pos[2][2];
    EKF[mdl_idx].P[2][3] = EKF[mdl_idx].P[3][2] = Pnew_pos[2][3];
    EKF[mdl_idx].P[2][4] = EKF[mdl_idx].P[4][2] = Pnew_pos[2][4];
    EKF[mdl_idx].P[3][3] = Pnew_pos[3][3];
    EKF[mdl_idx].P[3][4] = EKF[mdl_idx].P[4][3] = Pnew_pos[3][4];
    EKF[mdl_idx].P[4][4] = Pnew_pos[4][4];

    const ftype min_var = 1e-6f;
    for (uint8_t i = 0; i < 5; i++) {
        EKF[mdl_idx].P[i][i] = fmaxF(EKF[mdl_idx].P[i][i], min_var);
    }

    return true;
}

void EKFGSF_yaw_5state::resetEKFGSF()
{
    memset(&GSF, 0, sizeof(GSF));
    vel_fuse_running = false;
    run_ekf_gsf = false;
    pos_meas_dt = 1.0f;
    last_pos_fuse_ms = 0;
    n_clips = 0;

    memset(&EKF, 0, sizeof(EKF));
    const ftype yaw_increment = M_2PI / (ftype)N_MODELS_EKFGSF;
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        // evenly space initial yaw estimates in the region between +-Pi
        EKF[mdl_idx].X[2] = -M_PI + (0.5f * yaw_increment) + ((ftype)mdl_idx * yaw_increment);

        // All filter models start with the same weight
        GSF.weights[mdl_idx] = 1.0f / (ftype)N_MODELS_EKFGSF;

        // Use half yaw interval for yaw uncertainty as that is the maximum that the best model can be away from truth
        GSF.yaw_variance = sq(0.5f * yaw_increment);
        EKF[mdl_idx].P[2][2] = GSF.yaw_variance;
    }
}

// Apply a body frame delta angle to the body to earth frame rotation matrix using a small angle approximation
Matrix3F EKFGSF_yaw_5state::updateRotMat(const Matrix3F &R, const Vector3F &g) const
{
    Matrix3F ret = R;
    ret[0][0] += R[0][1] * g[2] - R[0][2] * g[1];
    ret[0][1] += R[0][2] * g[0] - R[0][0] * g[2];
    ret[0][2] += R[0][0] * g[1] - R[0][1] * g[0];
    ret[1][0] += R[1][1] * g[2] - R[1][2] * g[1];
    ret[1][1] += R[1][2] * g[0] - R[1][0] * g[2];
    ret[1][2] += R[1][0] * g[1] - R[1][1] * g[0];
    ret[2][0] += R[2][1] * g[2] - R[2][2] * g[1];
    ret[2][1] += R[2][2] * g[0] - R[2][0] * g[2];
    ret[2][2] += R[2][0] * g[1] - R[2][1] * g[0];

    // Renormalise rows
    ftype rowLengthSq;
    for (uint8_t r = 0; r < 3; r++) {
        rowLengthSq = ret[r][0] * ret[r][0] + ret[r][1] * ret[r][1] + ret[r][2] * ret[r][2];
        if (is_positive(rowLengthSq)) {
            // Use linear approximation for inverse sqrt taking advantage of the row length being close to 1.0
            const ftype rowLengthInv = 1.5f - 0.5f * rowLengthSq;
            Vector3F &row = ret[r];
            row *= rowLengthInv;
        }
    }

    return ret;
}

// returns true if a yaw estimate is available.  yaw and its variance
// is returned, as well as the number of models which are *not* being
// used to snthesise the yaw.
bool EKFGSF_yaw_5state::getYawData(ftype &yaw, ftype &yawVariance, uint8_t *_n_clips) const
{
    if (!vel_fuse_running) {
        return false;
    }
    yaw = GSF.yaw;
    yawVariance = GSF.yaw_variance;
    if (_n_clips != nullptr) {
        *_n_clips = n_clips;
    }
    return true;
}

bool EKFGSF_yaw_5state::getVelInnovLength(ftype &velInnovLength) const
{
    if (!vel_fuse_running) {
        return false;
    }
    velInnovLength = 0.0f;
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        velInnovLength += GSF.weights[mdl_idx] * sqrtF((sq(EKF[mdl_idx].innov_vel[0]) + sq(EKF[mdl_idx].innov_vel[1])));
    }
    return true;
}

void EKFGSF_yaw_5state::setGyroBias(Vector3f &gyroBias)
{
    for (uint8_t mdl_idx = 0; mdl_idx < N_MODELS_EKFGSF; mdl_idx++) {
        AHRS[mdl_idx].gyro_bias = gyroBias.toftype();
    }
}
