"""
simulation.py — синтетический симулятор 9-осевого IMU с ground truth.

Методика (агент methodology-advisor):
  * ground truth строится на мелкой сетке 1 кГц из гладких профилей углов Эйлера;
  * истинная угловая скорость тела вычисляется из кватернионной разности
    (кинематически согласована с ориентацией);
  * сенсоры сэмплируются на неравномерной сетке ~250 Гц с джиттером и пропусками;
  * все фильтры получают ОДНИ И ТЕ ЖЕ данные (одинаковый seed, шумы, помехи).

Сегменты сценария 'test' (180 с):
  S0   0–20   покой (калибровка bias)
  S1  20–40   медленные наклоны
  S2  40–55   быстрые повороты (до ~300 °/с)
  S3  55–70   рывки + динамическое линейное ускорение (до 6 м/с²)
  S4  70–85   вибрация (угловая 27 Гц + линейная 37/23 Гц)
  S5  85–105  магнитная помеха (вектор до 0.6·|B| + изменение амплитуды)
  S6 105–120  пропадание магнитометра (NaN)
  S7 120–150  комбинированное движение + скачок bias гироскопа на 130 с
  S8 150–180  замедление и возврат в покой

Явные допущения:
  * линейное ускорение задаётся напрямую в земной СК (позиция не интегрируется —
    для оценки ориентации важна только удельная сила);
  * hard/soft-iron носителя считаются откалиброванными; помеха — внешняя аддитивная;
  * при децимации (проверка устойчивости к dt) антиалиасинг не применяется —
    это ужесточает тест одинаково для всех фильтров.
"""
import numpy as np
from imu_math import euler_to_quat, quat_mult_batch, quat_to_euler_batch

DEG = np.pi / 180.0
G0 = 9.81
# ENU: магнитное поле севернее и вниз (наклонение 60°), нормированная амплитуда 1.0
MAG_DIP = 60.0 * DEG
MAG_REF_E = np.array([0.0, np.cos(MAG_DIP), -np.sin(MAG_DIP)])

SEGMENTS = {
    'static0':   (2.0, 20.0),
    'slow':      (20.0, 40.0),
    'fast':      (40.0, 55.0),
    'dyn_acc':   (55.0, 70.0),
    'vibration': (70.0, 85.0),
    'mag_dist':  (85.0, 105.0),
    'mag_out':   (105.0, 120.0),
    'combined':  (120.0, 150.0),
    'recovery':  (150.0, 178.0),
}
T_END = 180.0


def _win(t, t0, t1, rise=1.0):
    """Гладкое окно [t0, t1] с косинусными фронтами длительностью rise."""
    up = np.clip((t - t0) / rise, 0.0, 1.0)
    dn = np.clip((t1 - t) / rise, 0.0, 1.0)
    s = lambda x: 0.5 - 0.5 * np.cos(np.pi * x)
    return s(up) * s(dn)


def _euler_profiles(t, v):
    """Профили roll/pitch/yaw (рад) на мелкой сетке. v — словарь вариаций."""
    fs = v['freq_scale']; As = v['amp_scale']
    roll = np.zeros_like(t); pitch = np.zeros_like(t); yaw = np.zeros_like(t)

    w1 = _win(t, 20, 40, 2.0)     # slow
    roll += As*20*DEG*np.sin(2*np.pi*0.10*fs*(t-20)) * w1
    pitch += As*15*DEG*np.sin(2*np.pi*0.13*fs*(t-20) + 1.0) * w1
    yaw += As*10*DEG*np.sin(2*np.pi*0.05*fs*(t-20)) * w1

    w2 = _win(t, 40, 55, 1.0)     # fast
    yaw += As*120*DEG*np.sin(2*np.pi*0.25*fs*(t-40)) * w2
    roll += As*30*DEG*np.sin(2*np.pi*0.70*fs*(t-40)) * w2
    pitch += As*20*DEG*np.sin(2*np.pi*0.45*fs*(t-40) + 0.5) * w2

    w3 = _win(t, 55, 70, 0.8)     # jerky
    pitch += As*25*DEG*np.tanh(3*np.sin(2*np.pi*0.30*fs*(t-55))) * w3
    roll += As*15*DEG*np.tanh(4*np.sin(2*np.pi*0.45*fs*(t-55)+0.7)) * w3

    w4 = _win(t, 70, 85, 1.0)     # vibration: угловая компонента
    roll += 1.2*DEG*np.sin(2*np.pi*27*(t-70)) * w4
    pitch += 0.8*DEG*np.sin(2*np.pi*31*(t-70)+0.3) * w4
    yaw += As*15*DEG*np.sin(2*np.pi*0.08*fs*(t-70)) * w4

    w5 = _win(t, 85, 105, 1.5)    # во время магнитной помехи — умеренное движение
    roll += As*15*DEG*np.sin(2*np.pi*0.20*fs*(t-85)) * w5
    yaw += As*40*DEG*np.sin(2*np.pi*0.10*fs*(t-85)) * w5

    w6 = _win(t, 105, 120, 1.5)   # при пропадании магнитометра
    yaw += As*50*DEG*np.sin(2*np.pi*0.12*fs*(t-105)) * w6
    pitch += As*10*DEG*np.sin(2*np.pi*0.25*fs*(t-105)) * w6

    w7 = _win(t, 120, 150, 1.5)   # combined
    roll += As*25*DEG*np.sin(2*np.pi*0.30*fs*(t-120)) * w7
    pitch += As*20*DEG*np.sin(2*np.pi*0.22*fs*(t-120)+0.9) * w7
    yaw += As*60*DEG*np.sin(2*np.pi*0.15*fs*(t-120)) * w7

    w8 = _win(t, 150, 170, 2.0)   # recovery -> покой
    roll += As*8*DEG*np.sin(2*np.pi*0.10*fs*(t-150)) * w8
    return roll, pitch, yaw


def _quat_series(roll, pitch, yaw):
    Q = np.empty((len(roll), 4))
    for i in range(len(roll)):
        Q[i] = euler_to_quat(roll[i], pitch[i], yaw[i])
    # непрерывность знака
    for i in range(1, len(Q)):
        if np.dot(Q[i], Q[i-1]) < 0:
            Q[i] = -Q[i]
    return Q


def _body_rates(Q, dt_fine):
    """ω_body из кватернионной разности: q_{i+1} = q_i ⊗ exp(0.5 ω dt)."""
    Qc = Q[:-1].copy()
    Qc[:, 1:] *= -1.0                       # сопряжённые
    D = quat_mult_batch(Qc, Q[1:])          # q_i^{-1} ⊗ q_{i+1}
    D[D[:, 0] < 0] *= -1.0
    vecn = np.linalg.norm(D[:, 1:], axis=1)
    ang = 2.0 * np.arctan2(vecn, D[:, 0])
    scale = np.where(vecn > 1e-12, ang / np.maximum(vecn, 1e-12), 2.0)
    W = D[:, 1:] * scale[:, None] / dt_fine
    W = np.vstack([W, W[-1]])
    return W


def _dyn_accel_earth(t, rng, v):
    """Линейное динамическое ускорение в земной СК."""
    a = np.zeros((len(t), 3))
    # S3: импульсы каждые ~2 с, амплитуда до 6 м/с², случайные направления
    w3 = _win(t, 55, 70, 0.5)
    for k, t0 in enumerate(np.arange(56.0, 69.0, 2.0)):
        d = rng.standard_normal(3); d /= np.linalg.norm(d)
        pulse = 6.0 * np.exp(-((t - t0) / 0.18)**2)
        a += np.outer(pulse, d) * w3[:, None] * v['amp_scale']
    # S4: линейная вибрация
    w4 = _win(t, 70, 85, 0.8)
    a[:, 0] += 3.0*np.sin(2*np.pi*37*t) * w4
    a[:, 1] += 2.0*np.sin(2*np.pi*23*t + 0.4) * w4
    # S7: редкие умеренные толчки
    w7 = _win(t, 120, 150, 1.0)
    for t0 in (124.0, 133.0, 141.0, 147.0):
        d = rng.standard_normal(3); d /= np.linalg.norm(d)
        a += np.outer(3.0*np.exp(-((t - t0)/0.3)**2), d) * w7[:, None]
    return a


def _mag_disturbance_earth(t, v):
    """Аддитивная магнитная помеха в земной СК + множитель амплитуды."""
    w5 = _win(t, 85, 105, 2.0)
    dvec = np.array(v['mag_dist_vec'])      # доли |B|
    dist = np.outer(w5, dvec)
    amp = 1.0 + 0.25 * w5                   # изменение амплитуды поля
    return dist, amp


def generate_scenario(seed=42, variant='test', fs_nominal=250.0):
    """Возвращает словарь с данными сенсоров и ground truth."""
    rng = np.random.default_rng(seed)
    if variant == 'test':
        v = dict(freq_scale=1.0, amp_scale=1.0,
                 bias0=np.array([1.2, -0.8, 0.6]) * DEG,          # рад/с
                 bias_step=np.array([0.0, 0.4, 0.0]) * DEG,       # скачок на 130 с
                 mag_dist_vec=[0.45, -0.35, 0.25])
    elif variant == 'tune':
        v = dict(freq_scale=1.3, amp_scale=0.8,
                 bias0=np.array([-0.9, 1.1, -0.4]) * DEG,
                 bias_step=np.array([0.3, 0.0, -0.2]) * DEG,
                 mag_dist_vec=[-0.30, 0.40, -0.35])
    else:
        raise ValueError(variant)

    # ---- ground truth на мелкой сетке 1 кГц
    dt_f = 1e-3
    tf = np.arange(0.0, T_END + dt_f, dt_f)
    roll, pitch, yaw = _euler_profiles(tf, v)
    Qf = _quat_series(roll, pitch, yaw)
    Wf = _body_rates(Qf, dt_f)
    Af_dyn_e = _dyn_accel_earth(tf, rng, v)
    Mdist_e, Mamp = _mag_disturbance_earth(tf, v)

    # ---- неравномерная сетка сенсоров + пропуски данных
    dts = []
    tcur = 0.0
    gaps = [(62.0, 62.25), (95.0, 95.20), (140.0, 140.30)]
    base_dt = 1.0 / fs_nominal
    while tcur < T_END - base_dt:
        dt = base_dt * (1.0 + 0.10 * np.sin(2*np.pi*0.03*tcur)
                        + rng.uniform(-0.08, 0.08))
        tnew = tcur + dt
        for g0, g1 in gaps:
            if tcur < g0 <= tnew:
                tnew = g1                    # пропуск блока данных
        dts.append(tnew - tcur)
        tcur = tnew
    t = np.cumsum(dts)
    dt_arr = np.asarray(dts)
    N = len(t)

    # ---- интерполяция истины на сенсорные моменты
    idx = np.clip(np.searchsorted(tf, t) - 1, 0, len(tf) - 2)
    frac = ((t - tf[idx]) / dt_f)[:, None]
    Qt = Qf[idx] * (1 - frac) + Qf[idx + 1] * frac
    # знаковая согласованность при lerp
    flip = np.sum(Qf[idx] * Qf[idx + 1], axis=1) < 0
    Qt[flip] = (Qf[idx] * (1 - frac) - Qf[idx + 1] * frac)[flip]
    Qt /= np.linalg.norm(Qt, axis=1, keepdims=True)
    Wt = Wf[idx] * (1 - frac) + Wf[idx + 1] * frac
    Adyn_e = Af_dyn_e[idx] * (1 - frac) + Af_dyn_e[idx + 1] * frac
    Mdist = Mdist_e[idx] * (1 - frac) + Mdist_e[idx + 1] * frac
    MampS = (Mamp[idx] * (1 - frac[:, 0]) + Mamp[idx + 1] * frac[:, 0])

    # ---- модель гироскопа: bias (константа + случайное блуждание + скачок) + шум
    sigma_g = 0.15 * DEG                     # рад/с на отсчёт
    sigma_bw = 0.02 * DEG                    # рад/с/√с — случайное блуждание bias
    bias = np.zeros((N, 3))
    b = v['bias0'].copy()
    for i in range(N):
        b = b + rng.standard_normal(3) * sigma_bw * np.sqrt(dt_arr[i])
        bias[i] = b
        if t[i] >= 130.0 and (i == 0 or t[i-1] < 130.0):
            b = b + v['bias_step']
    gyr = Wt + bias + rng.standard_normal((N, 3)) * sigma_g

    # ---- модель акселерометра: R^T(g·ez + a_dyn) + шум (+вибрационный шум)
    sigma_a = 0.04
    w4s = _win(t, 70, 85, 0.8)
    acc = np.empty((N, 3))
    g_e = np.array([0.0, 0.0, G0])
    for i in range(N):
        w, x, y, z = Qt[i]
        R = np.array([
            [1-2*(y*y+z*z), 2*(x*y-w*z),   2*(x*z+w*y)],
            [2*(x*y+w*z),   1-2*(x*x+z*z), 2*(y*z-w*x)],
            [2*(x*z-w*y),   2*(y*z+w*x),   1-2*(x*x+y*y)]])
        acc[i] = R.T @ (g_e + Adyn_e[i])
    acc += rng.standard_normal((N, 3)) * (sigma_a + 0.5 * w4s[:, None])

    # ---- модель магнитометра: R^T((m_ref + помеха)·amp) + шум; NaN при пропадании
    sigma_m = 0.01
    mag = np.empty((N, 3))
    for i in range(N):
        w, x, y, z = Qt[i]
        R = np.array([
            [1-2*(y*y+z*z), 2*(x*y-w*z),   2*(x*z+w*y)],
            [2*(x*y+w*z),   1-2*(x*x+z*z), 2*(y*z-w*x)],
            [2*(x*z-w*y),   2*(y*z+w*x),   1-2*(x*x+y*y)]])
        m_e = (MAG_REF_E + Mdist[i]) * MampS[i]
        mag[i] = R.T @ m_e
    mag += rng.standard_normal((N, 3)) * sigma_m
    out_mask = (t >= 105.0) & (t < 120.0)
    mag[out_mask] = np.nan

    return dict(
        t=t, dt=dt_arr, gyr=gyr, acc=acc, mag=mag,
        q_true=Qt, euler_true=quat_to_euler_batch(Qt),
        gyro_bias_true=bias, segments=SEGMENTS,
        params=dict(g=G0, mag_ref_e=MAG_REF_E, sigma_g=sigma_g,
                    sigma_a=sigma_a, sigma_m=sigma_m, fs=fs_nominal),
    )


def decimate_scenario(data, factor):
    """Прореживание данных (проверка устойчивости к большим dt)."""
    sl = slice(None, None, factor)
    out = dict(data)
    for k in ('t', 'gyr', 'acc', 'mag', 'q_true', 'euler_true', 'gyro_bias_true'):
        out[k] = data[k][sl]
    tt = out['t']
    out['dt'] = np.diff(np.concatenate([[max(tt[0] - (tt[1]-tt[0]), 0.0)], tt]))
    return out


if __name__ == '__main__':
    d = generate_scenario()
    print(f"samples={len(d['t'])}, T={d['t'][-1]:.1f}s, "
          f"mean dt={d['dt'].mean()*1e3:.2f}ms, max dt={d['dt'].max()*1e3:.1f}ms")
    print(f"mag NaN: {np.isnan(d['mag'][:,0]).sum()} samples")
    print(f"|gyr| max={np.linalg.norm(d['gyr'],axis=1).max()/DEG:.0f} deg/s")
