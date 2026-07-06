"""Санити-проверки конвенций и сходимости (агент reason, фаза Induction-0)."""
import numpy as np
from math import pi, sin, cos, sqrt
from imu_math import (quat_mult, quat_normalize, quat_from_rotvec, quat_rotate,
                      quat_rotate_inv, gravity_body, quat_to_euler,
                      euler_to_quat, quat_angle)
from filters_baseline import (ComplementaryFilter, MahonyFilter,
                              MadgwickFilter, MEKF, init_quat_from_accel_mag)
from new_filter import ArbiterFilter

DEG = pi/180
G = 9.81
MAG_E = (0.0, 0.5, -0.8660254)
rng = np.random.default_rng(0)
fails = []


def check(name, cond, info=""):
    print(f"{'PASS' if cond else 'FAIL'}  {name} {info}")
    if not cond:
        fails.append(name)


# 1. roundtrip euler <-> quat
ok = True
for _ in range(200):
    r, p, y = rng.uniform(-pi, pi), rng.uniform(-pi/2*0.99, pi/2*0.99), rng.uniform(-pi, pi)
    q = euler_to_quat(r, p, y)
    r2, p2, y2 = quat_to_euler(q)
    if max(abs(r-r2), abs(p-p2), abs(y-y2)) > 1e-9:
        ok = False
check("euler<->quat roundtrip", ok)

# 2. quat_rotate соответствует R(q), gravity_body = R^T e_z
q = quat_normalize(tuple(rng.standard_normal(4)))
v = tuple(rng.standard_normal(3))
w, x, yq, z = q
R = np.array([[1-2*(yq*yq+z*z), 2*(x*yq-w*z), 2*(x*z+w*yq)],
              [2*(x*yq+w*z), 1-2*(x*x+z*z), 2*(yq*z-w*x)],
              [2*(x*z-w*yq), 2*(yq*z+w*x), 1-2*(x*x+yq*yq)]])
check("quat_rotate == R @ v", np.allclose(quat_rotate(q, v), R @ v, atol=1e-12))
check("gravity_body == R^T e_z", np.allclose(gravity_body(q), R.T @ [0, 0, 1], atol=1e-12))

# 3. кинематика: q <- q ⊗ exp(w dt) интегрирует вращение вокруг body-z
qq = (1.0, 0.0, 0.0, 0.0)
for _ in range(1000):
    qq = quat_normalize(quat_mult(qq, quat_from_rotvec(0, 0, 90*DEG*1e-3)))
_, _, yaw = quat_to_euler(qq)
check("gyro z-integration -> yaw 90deg", abs(yaw - 90*DEG) < 1e-6, f"yaw={yaw/DEG:.3f}")

# 4. init_quat_from_accel_mag восстанавливает известную ориентацию
q_true = euler_to_quat(20*DEG, -35*DEG, 140*DEG)
acc = quat_rotate_inv(q_true, (0, 0, G))
mag = quat_rotate_inv(q_true, MAG_E)
q0 = init_quat_from_accel_mag(acc, mag)
check("init from acc+mag", quat_angle(q0, q_true) < 0.2*DEG,
      f"err={quat_angle(q0, q_true)/DEG:.3f}deg")

# 5. статика с шумом: все фильтры <1.5deg через 20 с (старт с ошибкой 30deg yaw)
def run_static(f, T=60.0, dt=1/250):
    q_true = euler_to_quat(10*DEG, 5*DEG, 60*DEG)
    acc0 = quat_rotate_inv(q_true, (0, 0, G))
    mag0 = quat_rotate_inv(q_true, MAG_E)
    # намеренно портим инициализацию: первый отсчёт магнитометра повёрнут
    # на 30° вокруг земной вертикали -> начальная ошибка yaw = 30°
    mag_bad = quat_rotate_inv(quat_mult(quat_from_rotvec(0, 0, 30*DEG), q_true),
                              MAG_E)
    f.update((0.0, 0.0, 0.0), acc0, mag_bad, dt)
    n = int(T/dt)
    for i in range(n):
        gyr = tuple(rng.standard_normal(3)*0.1*DEG)
        acc = tuple(np.array(acc0) + rng.standard_normal(3)*0.03)
        mag = tuple(np.array(mag0) + rng.standard_normal(3)*0.01)
        f.update(gyr, acc, mag, dt)
    return quat_angle(f.q, q_true)/DEG

for F in (ComplementaryFilter(), MahonyFilter(), MadgwickFilter(),
          MEKF(), ArbiterFilter()):
    e = run_static(F)
    check(f"static+30deg-init-error converge [{F.name}]", e < 1.5, f"err={e:.2f}deg")

# 6. чистое вращение с идеальными датчиками: ARBITER не портит гироскоп
q_t = (1.0, 0.0, 0.0, 0.0)
f = ArbiterFilter()
dt = 1/250
err_max = 0.0
for i in range(int(30/dt)):
    t = i*dt
    w_b = (40*DEG*sin(2*pi*0.3*t), 30*DEG*cos(2*pi*0.2*t), 60*DEG*sin(2*pi*0.15*t))
    q_t = quat_normalize(quat_mult(q_t, quat_from_rotvec(w_b[0]*dt, w_b[1]*dt, w_b[2]*dt)))
    acc = quat_rotate_inv(q_t, (0, 0, G))
    mag = quat_rotate_inv(q_t, MAG_E)
    f.update(w_b, acc, mag, dt)
    if t > 2.0:
        err_max = max(err_max, quat_angle(f.q, q_t)/DEG)
check("ARBITER pure rotation err<1deg", err_max < 1.0, f"max={err_max:.3f}deg")

# 7. bias-обучение в статике: ARBITER находит bias 1.2 deg/s
f = ArbiterFilter()
q_true = euler_to_quat(0, 0, 0)
bias = np.array([1.2, -0.8, 0.6])*DEG
for i in range(int(20/dt)):
    gyr = tuple(bias + rng.standard_normal(3)*0.15*DEG)
    acc = tuple(np.array([0, 0, G]) + rng.standard_normal(3)*0.04)
    mag = tuple(np.array(MAG_E) + rng.standard_normal(3)*0.01)
    f.update(gyr, acc, mag, dt)
b_est = np.array([f.bx, f.by, f.bz])/DEG
check("ARBITER static bias est", np.max(np.abs(b_est - [1.2, -0.8, 0.6])) < 0.15,
      f"b_est={np.round(b_est,3)}")

# 8. симулятор: интегрирование истинного ω воспроизводит q_true
from simulation import generate_scenario
d = generate_scenario(seed=7)
gyr_clean = d['gyr'] - d['gyro_bias_true']
# уберём шум нельзя — проверяем на другом: интеграл ω_true из _body_rates
from simulation import _euler_profiles, _quat_series, _body_rates
tf = np.arange(0, 30, 1e-3)
v = dict(freq_scale=1.0, amp_scale=1.0)
r_, p_, y_ = _euler_profiles(tf, v)
Qf = _quat_series(r_, p_, y_)
Wf = _body_rates(Qf, 1e-3)
qi = tuple(Qf[0])
emax = 0.0
for i in range(len(tf)-1):
    qi = quat_normalize(quat_mult(qi, quat_from_rotvec(*(Wf[i]*1e-3))))
    emax = max(emax, quat_angle(qi, tuple(Qf[i+1])))
check("simulator omega consistency", emax/DEG < 0.5, f"max={emax/DEG:.4f}deg")

print("\n" + ("ALL OK" if not fails else f"FAILURES: {fails}"))
