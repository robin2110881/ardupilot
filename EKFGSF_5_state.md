# Goal
The 5 state gsf is design to provide a yaw source in the case of indoor compass less, gps denied environment where only cheap positionning system are used like uwb beacon

# How to use
SITL : 
1. launch ardupilot with gazebo
2. Add the rangefinder and optical flow sensor
3. Change EKF source params : change EK3_SRC1_ : POSXY to 8, POSZ to 2, VELXY to 5, VELZ to 0, YAW to 8.
4. Change to use EK3 5 state : EK3_SRC_OPTIONS to 16 (bit 4 to 1) 
5. Change VISO_OPTION to 1 to bypass the attitude reading from the vision_position_estimate messages

You can now lanch a program like this one to send horizontal position from gazebo to simulate a cheap position sensor, and to send startup yaw alignement (by pressing 'r' then enter), also you will need to set ekf origin by pressing 'o' : 


import os
import time
import traceback
import math
import argparse

from gz.transport13 import Node
from gz.msgs10.pose_v_pb2 import Pose_V
from pymavlink import mavutil

# --------------------------------
# Paramètres pouvant être modifiés
# --------------------------------

# Position HOME + Origine EKF
LATITUDE_HOME = -31.9508512 # Latitude HOME en deg
LONGITUDE_HOME = 115.863278 # Longitude HOME en deg
ALTITUDE_HOME = 10 # Altitude HOME en m

# Matrice de covariance (upper-triangular 6x6, must stay finite values for ArduPilot)
POSE_COV = [0.0 for _ in range(21)]
# Position sigma = 0.1m
POSE_COV[0] = 0.1**2   # x
POSE_COV[6] = 0.1**2   # y

POSE_TOPIC = "/world/iris_runway/dynamic_pose/info" # Topic Gazebo des positions

ODOMETRY_RATE = 12 # 12 Hz

SOURCE_COMPONENT = mavutil.mavlink.MAV_COMP_ID_VISUAL_INERTIAL_ODOMETRY

def address_conn_from_id(drone_id: int):
    return f"127.0.0.1:{14552 + drone_id*10}"

# --------------
# Initialisation
# --------------

# Eviter de recalculer à chaque fois
LATITUDE_HOME_E7 = int(LATITUDE_HOME * 1e7)
LONGITUDE_HOME_E7 = int(LONGITUDE_HOME * 1e7)
ALTITUDE_HOME_E3 = int(ALTITUDE_HOME * 1e3)

t_init = time.time() # Temps de boot

# On parse les arguments
parser = argparse.ArgumentParser(description="VisOdom for Gazebo")
parser.add_argument(
    "--world",
    required=True,
    type=str,
    help="Monde utilisé"
)
args = parser.parse_args()

# Nom du monde
WORLD = str(args.world)

# Version MAVLink + dialecte utilisé
os.environ['MAVLINK20'] = '1' # On met MAVLink 2
mavutil.set_dialect("ardupilotmega")

# Liste des ids (!= sysid) des drones disponibles
DRONE_IDS: list[int] = []
DRONE_IDS = [0]

# COnnexion MAVLink
def conn_with_drone(drone_id: int) -> mavutil.mavfile:
    """Renvoie une connexion MAVLink

    Args:
        drone_id (int): ID du drone (0 à 9)

    Returns:
        mavutil.mavfile: Connexion MAVLink
    """
    sysid = 11
    master = mavutil.mavlink_connection(
        device=address_conn_from_id(drone_id=drone_id),
        source_system=sysid,
        source_component=SOURCE_COMPONENT
    )
    master.first_byte = False # On set déjà MAVLink 2 au début
    master.target_system = sysid
    master.target_component = 1
    return master

# nom -> connexion
CONN_DICT: dict[str, mavutil.mavfile] = {f"iris_with_gimbal{k}": conn_with_drone(k) for k in DRONE_IDS}

# Classe pour gérer la position (évite d'arrondir dès qu'on reçoit la pos de Gazebo par ex)
class Position:
    def __init__(self):
        self._x = 0
        self._y = 0
        self._z = 0
        self.yaw = 0.0
        self.time_received = 0

    @property
    def x(self):
        return round(self._x, 2)

    @property
    def y(self):
        return round(self._y, 2)

    @property
    def z(self):
        return round(self._z, 2)

    @x.setter
    def x(self, other):
        self._x = other

    @y.setter
    def y(self, other):
        self._y = other

    @z.setter
    def z(self, other):
        self._z = other

    def is_valid(self) -> bool:
        return bool(time.time() - self.time_received <= 0.3)

# dict nom du drone -> position
dict_pos: dict[str, Position] = {}
for name in CONN_DICT:
    dict_pos[name] = Position()

import numpy as np
def pose_v_callback(msg) -> None:
    """Callback appelé à chaque message de type Pose_V reçu."""
    for pose in msg.pose:
        if pose.name in dict_pos:
            dict_pos[pose.name].x = pose.position.x
            dict_pos[pose.name].y = pose.position.y
            dict_pos[pose.name].z = pose.position.z
            # Convert Gazebo quaternion to yaw (rad)
            qx = pose.orientation.x
            qy = pose.orientation.y
            qz = pose.orientation.z
            qw = pose.orientation.w
            siny_cosp = 2.0 * (qw * qz + qx * qy)
            cosy_cosp = 1.0 - 2.0 * (qy * qy + qz * qz)
            # compute yaw and add small noise
            dict_pos[pose.name].yaw =  - math.atan2(siny_cosp, cosy_cosp) + math.pi/2 * (1+np.random.normal(0, 0.12)) # Gazebo a une orientation différente de celle d'ArduPilot, on aligne les axes
            dict_pos[pose.name].time_received = time.time()


# Connexion au topic
node = Node()
subscribed = node.subscribe(Pose_V, POSE_TOPIC, pose_v_callback)
if not subscribed:
    raise RuntimeError(f"Failed to subscribe to topic: {POSE_TOPIC}")

def send_odometry(conn_name: str):
    """Envoie un ODOMETRY

    Args:
        conn_name (str): Nom connexionHIGH_COV
    """
    if not dict_pos[conn_name].is_valid():
        return

    conn = CONN_DICT[conn_name]
    z_ned = -dict_pos[conn_name].z
    yaw = dict_pos[conn_name].yaw
    if not (math.isfinite(z_ned) and math.isfinite(yaw)):
        return
    conn.mav.vision_position_estimate_send(
        int((time.time() - t_init) * 1e6),
        dict_pos[conn_name].y, dict_pos[conn_name].x, z_ned,
        0.0, 0.0, yaw,
        POSE_COV,
        0
    )


def set_home_and_ekf_origin(conn: mavutil.mavfile):
    """Set HOME + EKF Origin

    Args:
        conn (mavutil.mavfile): Connexion
    """
    try:
        conn.mav.command_int_send(
            conn.target_system, # target_system
            mavutil.mavlink.MAV_COMP_ID_AUTOPILOT1, # target_component
            mavutil.mavlink.MAV_FRAME_GLOBAL, # frame
            mavutil.mavlink.MAV_CMD_DO_SET_HOME, # command
            0, # current (unused)
            0, # autocontinue (unused)
            0, # param1 (Use Current)
            0.0, # param2 (Roll)
            0.0, # param3 (Pitch)
            0.0, # param4 (Yaw)
            LATITUDE_HOME_E7, # x
            LONGITUDE_HOME_E7, # y
            ALTITUDE_HOME # z
        )
        conn.mav.set_gps_global_origin_send(
            conn.target_system, # target_system
            LATITUDE_HOME_E7, # latitude
            LONGITUDE_HOME_E7, # longitude
            ALTITUDE_HOME_E3, # altitude
            int(time.time() * 1e6) # time_usec
        )
    except Exception:
        print(traceback.format_exc())

# Tâches périodiques
odometry_task = mavutil.periodic_event(ODOMETRY_RATE)
odometry_task.force()

# ---------
# Main loop
# ---------
import threading

print("VisOdom started")
import sys

# Thread-safe flag for 'r' keypress
import threading

trigger_yaw_reset = threading.Event()
trigger_ekf_origin_reset = threading.Event()

def input_thread():
    while True:
        try:
            user_input = input()
            if user_input.strip().lower() == 'r':
                trigger_yaw_reset.set()
            if user_input.strip().lower() == 'o':
                trigger_ekf_origin_reset.set()
        except EOFError:
            break
        except Exception:
            print(traceback.format_exc())

# Start input thread
input_thr = threading.Thread(target=input_thread, daemon=True)
input_thr.start()

try:
    while True:
        # Drain MAVLink RX without blocking so buffers don't grow unbounded.
        for conn in CONN_DICT.values():
            _ = conn.recv_msg()

        if odometry_task.trigger():
            for conn_name in CONN_DICT:
                send_odometry(conn_name)
        if trigger_ekf_origin_reset.is_set():
            for conn in CONN_DICT.values():
                print("EKF origin reset sent")
                set_home_and_ekf_origin(conn)
                trigger_ekf_origin_reset.clear()
        if trigger_yaw_reset.is_set():
            for conn_name in CONN_DICT:
                print("Yaw align sent (with random noise)")
                yaw = dict_pos[conn_name].yaw
                sigma_rad = math.radians(10.0)
                yaw_var = sigma_rad * sigma_rad
                for conn in CONN_DICT.values():
                    conn.mav.command_int_send(
                        conn.target_system,
                        conn.target_component,
                        mavutil.mavlink.MAV_FRAME_GLOBAL_RELATIVE_ALT,
                        43006,   # AP_MAV_CMD_EXTERNAL_YAW_ESTIMATE (local)
                        0,       # current
                        0,       # autocontinue
                        yaw, # param1
                        yaw_var, # param2
                        0.0,     # param3 unused
                        0.0,     # param4 unused
                        0, 0, 0  # x, y, z unused
                    )
            trigger_yaw_reset.clear()

        # Yield CPU so Gazebo callback thread can run reliably.
        time.sleep(0.01)
except Exception:
    print(traceback.format_exc())
except KeyboardInterrupt:
    pass
finally:
    print("\nStopping VidOdom...")
    for conn in CONN_DICT.values():
        try:
            conn.close()
        except Exception:
            pass
    print("VidOdom stopped")
