"""
imu_math.py — кватернионные и векторные утилиты, общие для всех фильтров.

Конвенции (фиксированы для всего проекта):
  * Земная система координат: ENU (x=восток, y=север, z=вверх).
  * Кватернион q = [w, x, y, z] (Hamilton), отображает body -> earth:
        v_e = R(q) @ v_b,   R(q) — матрица поворота из q.
  * Кинематика: dq/dt = 0.5 * q ⊗ [0, ω_body].
  * Углы Эйлера ZYX (yaw-pitch-roll), применяются одинаково к истине и оценкам.

Скалярные функции написаны на чистом Python (кортежи) — для 3-векторов и
кватернионов это быстрее numpy и один-в-один переносится на C.
"""
from math import sqrt, sin, cos, asin, acos, atan2
import numpy as np


# ---------------------------------------------------------------- quaternions
def quat_mult(a, b):
    """Hamilton product a ⊗ b."""
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (aw*bw - ax*bx - ay*by - az*bz,
            aw*bx + ax*bw + ay*bz - az*by,
            aw*by - ax*bz + ay*bw + az*bx,
            aw*bz + ax*by - ay*bx + az*bw)


def quat_conj(q):
    return (q[0], -q[1], -q[2], -q[3])


def quat_normalize(q):
    w, x, y, z = q
    n = sqrt(w*w + x*x + y*y + z*z)
    if n < 1e-12:
        return (1.0, 0.0, 0.0, 0.0)
    inv = 1.0 / n
    return (w*inv, x*inv, y*inv, z*inv)


def quat_from_rotvec(rx, ry, rz):
    """Точная экспонента: вектор поворота (рад) -> кватернион."""
    angle = sqrt(rx*rx + ry*ry + rz*rz)
    if angle < 1e-9:
        # ряд Тейлора, достаточно 1-го порядка
        return quat_normalize((1.0, 0.5*rx, 0.5*ry, 0.5*rz))
    half = 0.5 * angle
    s = sin(half) / angle
    return (cos(half), rx*s, ry*s, rz*s)


def quat_rotate(q, v):
    """v_e = R(q) v_b (поворот вектора из body в earth)."""
    w, x, y, z = q
    vx, vy, vz = v
    # t = 2 * (q_vec × v)
    tx = 2.0*(y*vz - z*vy)
    ty = 2.0*(z*vx - x*vz)
    tz = 2.0*(x*vy - y*vx)
    # v' = v + w*t + q_vec × t
    return (vx + w*tx + (y*tz - z*ty),
            vy + w*ty + (z*tx - x*tz),
            vz + w*tz + (x*ty - y*tx))


def quat_rotate_inv(q, v):
    """v_b = R(q)ᵀ v_e (поворот вектора из earth в body)."""
    return quat_rotate(quat_conj(q), v)


def gravity_body(q):
    """Ожидаемое направление гравитации в body-СК: R(q)ᵀ [0,0,1] (3-я строка R)."""
    w, x, y, z = q
    return (2.0*(x*z - w*y), 2.0*(y*z + w*x), 1.0 - 2.0*(x*x + y*y))


def quat_to_euler(q):
    """ZYX: roll (X), pitch (Y), yaw (Z), рад."""
    w, x, y, z = q
    sinp = 2.0*(w*y - z*x)
    sinp = max(-1.0, min(1.0, sinp))
    roll = atan2(2.0*(w*x + y*z), 1.0 - 2.0*(x*x + y*y))
    pitch = asin(sinp)
    yaw = atan2(2.0*(w*z + x*y), 1.0 - 2.0*(y*y + z*z))
    return roll, pitch, yaw


def euler_to_quat(roll, pitch, yaw):
    """Обратное преобразование, композиция qz(yaw) ⊗ qy(pitch) ⊗ qx(roll)."""
    cr, sr = cos(roll*0.5), sin(roll*0.5)
    cp, sp = cos(pitch*0.5), sin(pitch*0.5)
    cy, sy = cos(yaw*0.5), sin(yaw*0.5)
    w = cy*cp*cr + sy*sp*sr
    x = cy*cp*sr - sy*sp*cr
    y = cy*sp*cr + sy*cp*sr
    z = sy*cp*cr - cy*sp*sr
    return (w, x, y, z)


def quat_angle(qa, qb):
    """Геодезический угол между ориентациями, рад."""
    d = abs(qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2] + qa[3]*qb[3])
    d = min(1.0, d)
    return 2.0 * acos(d)


# ------------------------------------------------------------- batch (numpy)
def quat_to_euler_batch(Q):
    """Q: (N,4) -> (N,3) roll,pitch,yaw."""
    w, x, y, z = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
    sinp = np.clip(2.0*(w*y - z*x), -1.0, 1.0)
    roll = np.arctan2(2.0*(w*x + y*z), 1.0 - 2.0*(x*x + y*y))
    pitch = np.arcsin(sinp)
    yaw = np.arctan2(2.0*(w*z + x*y), 1.0 - 2.0*(y*y + z*z))
    return np.stack([roll, pitch, yaw], axis=1)


def quat_angle_batch(QA, QB):
    d = np.abs(np.sum(QA * QB, axis=1))
    return 2.0 * np.arccos(np.clip(d, -1.0, 1.0))


def quat_mult_batch(A, B):
    aw, ax, ay, az = A[:, 0], A[:, 1], A[:, 2], A[:, 3]
    bw, bx, by, bz = B[:, 0], B[:, 1], B[:, 2], B[:, 3]
    return np.stack([
        aw*bw - ax*bx - ay*by - az*bz,
        aw*bx + ax*bw + ay*bz - az*by,
        aw*by - ax*bz + ay*bw + az*bx,
        aw*bz + ax*by - ay*bx + az*bw], axis=1)


def wrap_angle(a):
    """Свёртка угла/массива углов в [-pi, pi]."""
    return (a + np.pi) % (2.0*np.pi) - np.pi
