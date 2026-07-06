"""
filters_baseline.py — эталонные алгоритмы для честного сравнения.

Реализованы по первоисточникам:
  * ComplementaryFilter — кватернионный комплементарный фильтр (константные
    коэффициенты, коррекция наклона по акселерометру и yaw по магнитометру);
  * MahonyFilter — explicit complementary filter (Mahony, Hamel, Pflimlin 2008)
    с PI-обратной связью и оценкой bias через интегральный член;
  * MadgwickFilter — градиентный фильтр (Madgwick 2010, MARG-версия,
    транскрипция эталонного кода x-io);
  * MEKF — мультипликативный расширенный фильтр Калмана, error-state 6
    (3 ошибки ориентации + 3 bias), без адаптации R (базовый вариант).

Единый API: f.update(gyr, acc, mag, dt) -> q (tuple w,x,y,z)
  gyr: (gx,gy,gz) рад/с; acc: (ax,ay,az) м/с²; mag: (mx,my,mz) или None; dt: с.
Все фильтры инициализируются одинаково: init_quat_from_accel_mag().
Магнитометр во всех фильтрах используется в azimuth-free форме (не влияет
склонение); NaN/None магнитометра корректно пропускается.
"""
from math import sqrt, atan2, acos, sin, cos
import numpy as np
from imu_math import (quat_mult, quat_normalize, quat_from_rotvec,
                      quat_rotate, gravity_body, quat_conj)


def init_quat_from_accel_mag(acc, mag=None):
    """Общая инициализация: наклон из акселерометра, yaw из магнитометра."""
    ax, ay, az = acc
    n = sqrt(ax*ax + ay*ay + az*az)
    if n < 1e-6:
        return (1.0, 0.0, 0.0, 0.0)
    ax, ay, az = ax/n, ay/n, az/n
    # кватернион, поворачивающий â (body) в e_z (earth): tilt-only
    cx, cy, cz = ay*1.0 - az*0.0, az*0.0 - ax*1.0, ax*0.0 - ay*0.0  # â × e_z
    s = sqrt(cx*cx + cy*cy + cz*cz)
    c = az  # â · e_z
    if s < 1e-9:
        q_t = (1.0, 0.0, 0.0, 0.0) if c > 0 else (0.0, 1.0, 0.0, 0.0)
    else:
        ang = atan2(s, c)
        q_t = quat_from_rotvec(cx/s*ang, cy/s*ang, cz/s*ang)
    if mag is None or any(m != m for m in mag):
        return q_t
    m_e = quat_rotate(q_t, mag)
    h = sqrt(m_e[0]*m_e[0] + m_e[1]*m_e[1])
    if h < 1e-6:
        return q_t
    # m_e имеет азимут psi; ошибка yaw = -psi, поэтому коррекция +psi
    psi = atan2(m_e[0], m_e[1])
    return quat_normalize(quat_mult(quat_from_rotvec(0.0, 0.0, psi), q_t))


def _mag_ok(mag):
    return mag is not None and mag[0] == mag[0] and mag[1] == mag[1] and mag[2] == mag[2]


# =========================================================== Complementary
class ComplementaryFilter:
    """Классический кватернионный комплементарный фильтр с константными
    постоянными времени tau_a (наклон) и tau_m (yaw)."""

    name = 'Complementary'

    def __init__(self, tau_a=2.0, tau_m=5.0):
        self.tau_a = tau_a
        self.tau_m = tau_m
        self.q = (1.0, 0.0, 0.0, 0.0)
        self._init = False

    def update(self, gyr, acc, mag, dt):
        if not self._init:
            self.q = init_quat_from_accel_mag(acc, mag if _mag_ok(mag) else None)
            self._init = True
            return self.q
        q = quat_mult(self.q, quat_from_rotvec(gyr[0]*dt, gyr[1]*dt, gyr[2]*dt))
        # --- коррекция наклона: доля k_a от рассогласования â и v̂
        ax, ay, az = acc
        n = sqrt(ax*ax + ay*ay + az*az)
        if n > 1e-6:
            ax, ay, az = ax/n, ay/n, az/n
            vx, vy, vz = gravity_body(q)
            ex, ey, ez = ay*vz - az*vy, az*vx - ax*vz, ax*vy - ay*vx  # â × v̂
            s = sqrt(ex*ex + ey*ey + ez*ez)
            c = ax*vx + ay*vy + az*vz
            if s > 1e-9:
                ang = atan2(s, c)
                k_a = dt / (self.tau_a + dt)
                f = k_a * ang / s
                q = quat_mult(q, quat_from_rotvec(ex*f, ey*f, ez*f))
        # --- коррекция yaw по магнитометру
        if _mag_ok(mag):
            m_e = quat_rotate(q, mag)
            h = sqrt(m_e[0]*m_e[0] + m_e[1]*m_e[1])
            if h > 1e-6:
                psi = atan2(m_e[0], m_e[1])   # ошибка yaw = -psi
                k_m = dt / (self.tau_m + dt)
                q = quat_mult(quat_from_rotvec(0.0, 0.0, k_m*psi), q)
        self.q = quat_normalize(q)
        return self.q


# ================================================================= Mahony
class MahonyFilter:
    """Explicit complementary filter (Mahony et al. 2008): PI-обратная связь
    по ошибке e = Σ v_meas × v_est, bias — через интегральный член Ki."""

    name = 'Mahony'

    def __init__(self, kp=1.0, ki=0.1):
        self.kp = kp
        self.ki = ki
        self.q = (1.0, 0.0, 0.0, 0.0)
        self.bx = self.by = self.bz = 0.0   # интегральный член (−bias)
        self._init = False

    def update(self, gyr, acc, mag, dt):
        if not self._init:
            self.q = init_quat_from_accel_mag(acc, mag if _mag_ok(mag) else None)
            self._init = True
            return self.q
        q = self.q
        ex = ey = ez = 0.0
        ax, ay, az = acc
        n = sqrt(ax*ax + ay*ay + az*az)
        if n > 1e-6:
            ax, ay, az = ax/n, ay/n, az/n
            vx, vy, vz = gravity_body(q)
            ex += ay*vz - az*vy
            ey += az*vx - ax*vz
            ez += ax*vy - ay*vx
        if _mag_ok(mag):
            mx, my, mz = mag
            nm = sqrt(mx*mx + my*my + mz*mz)
            if nm > 1e-6:
                mx, my, mz = mx/nm, my/nm, mz/nm
                # azimuth-free эталон: h = R m̂, b = [0, √(hx²+hy²), hz]
                hx, hy, hz = quat_rotate(q, (mx, my, mz))
                bref_y = sqrt(hx*hx + hy*hy)
                wx, wy, wz = quat_rotate(quat_conj(q), (0.0, bref_y, hz))
                ex += my*wz - mz*wy
                ey += mz*wx - mx*wz
                ez += mx*wy - my*wx
        if self.ki > 0.0:
            self.bx += self.ki * ex * dt
            self.by += self.ki * ey * dt
            self.bz += self.ki * ez * dt
        gx = gyr[0] + self.kp*ex + self.bx
        gy = gyr[1] + self.kp*ey + self.by
        gz = gyr[2] + self.kp*ez + self.bz
        self.q = quat_normalize(quat_mult(q, quat_from_rotvec(gx*dt, gy*dt, gz*dt)))
        return self.q


# =============================================================== Madgwick
class MadgwickFilter:
    """Градиентный фильтр Madgwick (2010), MARG-версия. Шаг — нормированный
    градиент функции ошибки, амплитуда шага — константа beta."""

    name = 'Madgwick'

    # Эталонный код x-io использует azimuth-free эталон b=[bx,0,bz] — его
    # земная СК NWU (север=+x). Для сопоставимости с ENU (север=+y) внутреннее
    # состояние держим в NWU и конвертируем на входе/выходе поворотом ±90° по z.
    _Q_NWU2ENU = (0.7071067811865476, 0.0, 0.0, 0.7071067811865476)
    _Q_ENU2NWU = (0.7071067811865476, 0.0, 0.0, -0.7071067811865476)

    def __init__(self, beta=0.05):
        self.beta = beta
        self._qi = (1.0, 0.0, 0.0, 0.0)      # внутренний кватернион (NWU)
        self.q = (1.0, 0.0, 0.0, 0.0)
        self._init = False

    def update(self, gyr, acc, mag, dt):
        if not self._init:
            q_enu = init_quat_from_accel_mag(acc, mag if _mag_ok(mag) else None)
            self._qi = quat_normalize(quat_mult(self._Q_ENU2NWU, q_enu))
            self.q = q_enu
            self._init = True
            return self.q
        q0, q1, q2, q3 = self._qi
        gx, gy, gz = gyr
        ax, ay, az = acc
        na = sqrt(ax*ax + ay*ay + az*az)
        use_mag = _mag_ok(mag)
        if na > 1e-6:
            ax, ay, az = ax/na, ay/na, az/na
        if use_mag:
            mx, my, mz = mag
            nm = sqrt(mx*mx + my*my + mz*mz)
            if nm > 1e-6:
                mx, my, mz = mx/nm, my/nm, mz/nm
            else:
                use_mag = False

        if na > 1e-6 and use_mag:
            # эталонная MARG-форма (x-io)
            _2q0mx = 2.0*q0*mx; _2q0my = 2.0*q0*my; _2q0mz = 2.0*q0*mz
            _2q1mx = 2.0*q1*mx
            _2q0 = 2.0*q0; _2q1 = 2.0*q1; _2q2 = 2.0*q2; _2q3 = 2.0*q3
            _2q0q2 = 2.0*q0*q2; _2q2q3 = 2.0*q2*q3
            q0q0 = q0*q0; q0q1 = q0*q1; q0q2 = q0*q2; q0q3 = q0*q3
            q1q1 = q1*q1; q1q2 = q1*q2; q1q3 = q1*q3
            q2q2 = q2*q2; q2q3 = q2*q3; q3q3 = q3*q3

            hx = mx*q0q0 - _2q0my*q3 + _2q0mz*q2 + mx*q1q1 + _2q1*my*q2 \
                + _2q1*mz*q3 - mx*q2q2 - mx*q3q3
            hy = _2q0mx*q3 + my*q0q0 - _2q0mz*q1 + _2q1mx*q2 - my*q1q1 \
                + my*q2q2 + _2q2*mz*q3 - my*q3q3
            _2bx = sqrt(hx*hx + hy*hy)
            _2bz = -_2q0mx*q2 + _2q0my*q1 + mz*q0q0 + _2q1mx*q3 - mz*q1q1 \
                + _2q2*my*q3 - mz*q2q2 + mz*q3q3
            _4bx = 2.0*_2bx; _4bz = 2.0*_2bz

            s0 = -_2q2*(2.0*q1q3 - _2q0q2 - ax) + _2q1*(2.0*q0q1 + _2q2q3 - ay) \
                - _2bz*q2*(_2bx*(0.5 - q2q2 - q3q3) + _2bz*(q1q3 - q0q2) - mx) \
                + (-_2bx*q3 + _2bz*q1)*(_2bx*(q1q2 - q0q3) + _2bz*(q0q1 + q2q3) - my) \
                + _2bx*q2*(_2bx*(q0q2 + q1q3) + _2bz*(0.5 - q1q1 - q2q2) - mz)
            s1 = _2q3*(2.0*q1q3 - _2q0q2 - ax) + _2q0*(2.0*q0q1 + _2q2q3 - ay) \
                - 4.0*q1*(1.0 - 2.0*q1q1 - 2.0*q2q2 - az) \
                + _2bz*q3*(_2bx*(0.5 - q2q2 - q3q3) + _2bz*(q1q3 - q0q2) - mx) \
                + (_2bx*q2 + _2bz*q0)*(_2bx*(q1q2 - q0q3) + _2bz*(q0q1 + q2q3) - my) \
                + (_2bx*q3 - _4bz*q1)*(_2bx*(q0q2 + q1q3) + _2bz*(0.5 - q1q1 - q2q2) - mz)
            s2 = -_2q0*(2.0*q1q3 - _2q0q2 - ax) + _2q3*(2.0*q0q1 + _2q2q3 - ay) \
                - 4.0*q2*(1.0 - 2.0*q1q1 - 2.0*q2q2 - az) \
                + (-_4bx*q2 - _2bz*q0)*(_2bx*(0.5 - q2q2 - q3q3) + _2bz*(q1q3 - q0q2) - mx) \
                + (_2bx*q1 + _2bz*q3)*(_2bx*(q1q2 - q0q3) + _2bz*(q0q1 + q2q3) - my) \
                + (_2bx*q0 - _4bz*q2)*(_2bx*(q0q2 + q1q3) + _2bz*(0.5 - q1q1 - q2q2) - mz)
            s3 = _2q1*(2.0*q1q3 - _2q0q2 - ax) + _2q2*(2.0*q0q1 + _2q2q3 - ay) \
                + (-_4bx*q3 + _2bz*q1)*(_2bx*(0.5 - q2q2 - q3q3) + _2bz*(q1q3 - q0q2) - mx) \
                + (-_2bx*q0 + _2bz*q2)*(_2bx*(q1q2 - q0q3) + _2bz*(q0q1 + q2q3) - my) \
                + _2bx*q1*(_2bx*(q0q2 + q1q3) + _2bz*(0.5 - q1q1 - q2q2) - mz)
        elif na > 1e-6:
            # IMU-версия (без магнитометра)
            _2q0 = 2.0*q0; _2q1 = 2.0*q1; _2q2 = 2.0*q2; _2q3 = 2.0*q3
            _4q0 = 4.0*q0; _4q1 = 4.0*q1; _4q2 = 4.0*q2
            _8q1 = 8.0*q1; _8q2 = 8.0*q2
            q0q0 = q0*q0; q1q1 = q1*q1; q2q2 = q2*q2; q3q3 = q3*q3
            s0 = _4q0*q2q2 + _2q2*ax + _4q0*q1q1 - _2q1*ay
            s1 = _4q1*q3q3 - _2q3*ax + 4.0*q0q0*q1 - _2q0*ay - _4q1 \
                + _8q1*q1q1 + _8q1*q2q2 + _4q1*az
            s2 = 4.0*q0q0*q2 + _2q0*ax + _4q2*q3q3 - _2q3*ay - _4q2 \
                + _8q2*q1q1 + _8q2*q2q2 + _4q2*az
            s3 = 4.0*q1q1*q3 - _2q1*ax + 4.0*q2q2*q3 - _2q2*ay
        else:
            s0 = s1 = s2 = s3 = 0.0

        qDot0 = 0.5*(-q1*gx - q2*gy - q3*gz)
        qDot1 = 0.5*(q0*gx + q2*gz - q3*gy)
        qDot2 = 0.5*(q0*gy - q1*gz + q3*gx)
        qDot3 = 0.5*(q0*gz + q1*gy - q2*gx)
        ns = sqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3)
        if ns > 1e-12:
            inv = self.beta / ns
            qDot0 -= s0*inv; qDot1 -= s1*inv; qDot2 -= s2*inv; qDot3 -= s3*inv
        self._qi = quat_normalize((q0 + qDot0*dt, q1 + qDot1*dt,
                                   q2 + qDot2*dt, q3 + qDot3*dt))
        self.q = quat_normalize(quat_mult(self._Q_NWU2ENU, self._qi))
        return self.q


# ==================================================================== MEKF
class MEKF:
    """Мультипликативный EKF, error-state x=[δθ(3), δb(3)], P 6×6.
    Базовый вариант без адаптации/гейтирования измерений — честный
    представитель EKF-класса для сравнения точности и стоимости."""

    name = 'MEKF'

    def __init__(self, sigma_g=0.003, sigma_bw=0.0005,
                 sigma_acc_dir=0.05, sigma_mag_dir=0.2,
                 p0_att=0.3, p0_bias=0.03):
        self.q = (1.0, 0.0, 0.0, 0.0)
        self.b = np.zeros(3)
        self.P = np.diag([p0_att**2]*3 + [p0_bias**2]*3)
        self.sg2 = sigma_g**2
        self.sbw2 = sigma_bw**2
        self.Ra = np.eye(3) * sigma_acc_dir**2
        self.Rm = np.eye(3) * sigma_mag_dir**2
        self._init = False
        self._I6 = np.eye(6)

    def _vec_update(self, v_meas, v_pred, R):
        H = np.zeros((3, 6))
        H[:, :3] = np.array([[0.0, -v_pred[2], v_pred[1]],
                             [v_pred[2], 0.0, -v_pred[0]],
                             [-v_pred[1], v_pred[0], 0.0]])
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.inv(S)
        dx = K @ (np.asarray(v_meas) - np.asarray(v_pred))
        self.P = (self._I6 - K @ H) @ self.P
        self.P = 0.5 * (self.P + self.P.T)
        return dx

    def update(self, gyr, acc, mag, dt):
        if not self._init:
            self.q = init_quat_from_accel_mag(acc, mag if _mag_ok(mag) else None)
            self._init = True
            return self.q
        w = np.asarray(gyr) - self.b
        self.q = quat_normalize(quat_mult(
            self.q, quat_from_rotvec(w[0]*dt, w[1]*dt, w[2]*dt)))
        # пропагация ковариации
        Wx = np.array([[0.0, -w[2], w[1]],
                       [w[2], 0.0, -w[0]],
                       [-w[1], w[0], 0.0]])
        F = self._I6.copy()
        F[:3, :3] -= Wx * dt
        F[:3, 3:] = -np.eye(3) * dt
        self.P = F @ self.P @ F.T
        self.P[:3, :3] += np.eye(3) * (self.sg2 * dt)
        self.P[3:, 3:] += np.eye(3) * (self.sbw2 * dt)

        dx = np.zeros(6)
        ax, ay, az = acc
        na = sqrt(ax*ax + ay*ay + az*az)
        if na > 1e-6:
            a_hat = (ax/na, ay/na, az/na)
            v_pred = gravity_body(self.q)
            dx = self._vec_update(a_hat, v_pred, self.Ra)
            self._inject(dx)
        if _mag_ok(mag):
            mx, my, mz = mag
            nm = sqrt(mx*mx + my*my + mz*mz)
            if nm > 1e-6:
                m_hat = (mx/nm, my/nm, mz/nm)
                hx, hy, hz = quat_rotate(self.q, m_hat)
                bref = (0.0, sqrt(hx*hx + hy*hy), hz)
                w_pred = quat_rotate(quat_conj(self.q), bref)
                dx = self._vec_update(m_hat, w_pred, self.Rm)
                self._inject(dx)
        return self.q

    def _inject(self, dx):
        self.q = quat_normalize(quat_mult(
            self.q, quat_from_rotvec(dx[0], dx[1], dx[2])))
        self.b = self.b + dx[3:]
