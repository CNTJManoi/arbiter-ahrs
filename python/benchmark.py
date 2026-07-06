"""
benchmark.py — честное сравнение алгоритмов на синтетических данных.

Протокол (агент methodology-advisor):
  1. ТЮНИНГ: сеточный поиск параметров КАЖДОГО фильтра на сценарии
     variant='tune' (seed=1: другие траектории, bias, вектор помехи).
     Критерий — RMSE полного угла ориентации (после 5 с переходного).
  2. ТЕСТ: один прогон на variant='test' (seed=42) с лучшими параметрами.
     Никакого перетюнинга на тестовых данных.
  3. Посегментные метрики, тайминг на шаг, устойчивость к dt (децимация
     до ~50 Гц и ~21 Гц с теми же параметрами).

Запуск:  python benchmark.py            (полный тюнинг + тест)
         python benchmark.py --no-tune  (параметры из results/tuned_params.json)
"""
import json
import sys
import time
import itertools
import numpy as np

from simulation import generate_scenario, decimate_scenario, SEGMENTS
from imu_math import quat_to_euler_batch, quat_angle_batch, wrap_angle
from filters_baseline import (ComplementaryFilter, MahonyFilter,
                              MadgwickFilter, MEKF)
from new_filter import ArbiterFilter

DEG = np.pi / 180.0
RESULTS = 'results'

FILTER_SPECS = {
    'Complementary': (ComplementaryFilter, dict(
        tau_a=[0.2, 0.5, 1.0, 2.0, 5.0], tau_m=[0.3, 0.5, 1.0, 2.0, 5.0, 10.0])),
    'Mahony': (MahonyFilter, dict(
        kp=[0.3, 0.5, 1.0, 2.0, 4.0, 8.0], ki=[0.0, 0.05, 0.1, 0.3, 1.0])),
    'Madgwick': (MadgwickFilter, dict(
        beta=[0.01, 0.02, 0.033, 0.05, 0.1, 0.2, 0.4, 0.8])),
    'MEKF': (MEKF, dict(
        sigma_acc_dir=[0.01, 0.02, 0.05, 0.15, 0.4],
        sigma_mag_dir=[0.01, 0.02, 0.05, 0.2, 0.5])),
    'ARBITER': (ArbiterFilter, dict(
        rho0=[0.075*DEG, 0.15*DEG, 0.3*DEG, 0.6*DEG],
        diff_tol=[0.5, 0.8, 1.2],
        cap_yaw=[5.0*DEG, 15.0*DEG],
        cap_tilt=[5.0*DEG, 10.0*DEG])),
}


def run_filter(f, data, record_diag=False):
    """Прогон фильтра по данным; возвращает (Q, per_step_us, diag|None)."""
    gyr, acc, mag, dt = data['gyr'], data['acc'], data['mag'], data['dt']
    N = len(dt)
    Q = np.empty((N, 4))
    diag = {k: np.zeros(N) for k in
            ('Ea', 'Em', 'bx', 'by', 'bz', 'Tt', 'Ty', 'ub',
             'static', 'gate_a', 'gate_m')} if record_diag else None
    t0 = time.perf_counter()
    for i in range(N):
        q = f.update((gyr[i, 0], gyr[i, 1], gyr[i, 2]),
                     (acc[i, 0], acc[i, 1], acc[i, 2]),
                     (mag[i, 0], mag[i, 1], mag[i, 2]),
                     dt[i])
        Q[i, 0], Q[i, 1], Q[i, 2], Q[i, 3] = q
        if record_diag:
            diag['Ea'][i] = f.Ea; diag['Em'][i] = f.Em
            diag['bx'][i] = f.bx; diag['by'][i] = f.by; diag['bz'][i] = f.bz
            diag['Tt'][i] = f.Tt; diag['Ty'][i] = f.Ty; diag['ub'][i] = f.ub
            diag['static'][i] = f.is_static
            diag['gate_a'][i] = f.gate_a; diag['gate_m'][i] = f.gate_m
    per_step_us = (time.perf_counter() - t0) / N * 1e6
    # выравнивание знака для метрик
    sign = np.where(np.sum(Q * data['q_true'], axis=1) < 0, -1.0, 1.0)
    Q *= sign[:, None]
    return Q, per_step_us, diag


def _tilt_err_deg(Q_est, Q_true):
    """Ошибка наклона: угол между оценкой и истиной направления гравитации."""
    def gb(Q):
        w, x, y, z = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
        return np.stack([2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y)], axis=1)
    d = np.clip(np.sum(gb(Q_est) * gb(Q_true), axis=1), -1.0, 1.0)
    return np.degrees(np.arccos(d))


def compute_metrics(Q_est, data, skip=5.0):
    t = data['t']
    m = t > skip
    ang = np.degrees(quat_angle_batch(Q_est, data['q_true']))
    tilt = _tilt_err_deg(Q_est, data['q_true'])
    eul = quat_to_euler_batch(Q_est)
    de = np.degrees(wrap_angle(eul - data['euler_true']))
    res = dict(
        rmse_roll=float(np.sqrt(np.mean(de[m, 0]**2))),
        rmse_pitch=float(np.sqrt(np.mean(de[m, 1]**2))),
        rmse_yaw=float(np.sqrt(np.mean(de[m, 2]**2))),
        rmse_ang=float(np.sqrt(np.mean(ang[m]**2))),
        mean_ang=float(np.mean(ang[m])),
        max_ang=float(np.max(ang[m])),
    )
    for name, (t0, t1) in SEGMENTS.items():
        ms = (t >= t0) & (t < t1)
        if ms.sum() < 10:
            continue
        res[f'seg_{name}_ang'] = float(np.sqrt(np.mean(ang[ms]**2)))
        res[f'seg_{name}_tilt'] = float(np.sqrt(np.mean(tilt[ms]**2)))
        res[f'seg_{name}_yaw'] = float(np.sqrt(np.mean(de[ms, 2]**2)))
    return res


def tune_all():
    data_tune = generate_scenario(seed=1, variant='tune')
    best = {}
    for name, (cls, grid) in FILTER_SPECS.items():
        keys = list(grid.keys())
        best_score, best_params = np.inf, {}
        for combo in itertools.product(*(grid[k] for k in keys)):
            params = dict(zip(keys, combo))
            Q, _, _ = run_filter(cls(**params), data_tune)
            score = compute_metrics(Q, data_tune)['rmse_ang']
            if score < best_score:
                best_score, best_params = score, params
        best[name] = dict(params=best_params, tune_rmse=best_score)
        print(f"  tuned {name:14s} rmse={best_score:6.2f} deg  {best_params}")
    with open(f'{RESULTS}/tuned_params.json', 'w') as fh:
        json.dump(best, fh, indent=2)
    return best


def main():
    if '--no-tune' in sys.argv:
        with open(f'{RESULTS}/tuned_params.json') as fh:
            best = json.load(fh)
    else:
        print("=== ТЮНИНГ (variant='tune', seed=1) ===")
        best = tune_all()

    print("\n=== ТЕСТ (variant='test', seed=42) ===")
    data = generate_scenario(seed=42, variant='test')
    all_metrics, trajs, timings = {}, {}, {}
    diag_saved = None
    for name, (cls, _) in FILTER_SPECS.items():
        params = best[name]['params']
        f = cls(**params)
        Q, us, diag = run_filter(f, data, record_diag=(name == 'ARBITER'))
        # тайминг: повторный прогон без записи диагностики для чистоты
        _, us2, _ = run_filter(cls(**params), data)
        us = min(us, us2)
        m = compute_metrics(Q, data)
        m['step_us'] = us
        all_metrics[name] = m
        trajs[name] = Q
        timings[name] = us
        if diag is not None:
            diag_saved = diag
        print(f"  {name:14s} rmse_ang={m['rmse_ang']:6.2f}  "
              f"tilt(dyn)={m['seg_dyn_acc_tilt']:5.2f}  "
              f"yaw(mag_dist)={m['seg_mag_dist_yaw']:6.2f}  "
              f"{us:6.1f} us/step")

    # --- устойчивость к dt (те же параметры, прореженные данные)
    print("\n=== УСТОЙЧИВОСТЬ К dt ===")
    for factor, label in ((5, '50Hz'), (12, '21Hz')):
        dec = decimate_scenario(data, factor)
        for name, (cls, _) in FILTER_SPECS.items():
            Q, _, _ = run_filter(cls(**best[name]['params']), dec)
            r = compute_metrics(Q, dec)['rmse_ang']
            all_metrics[name][f'rmse_ang_{label}'] = r
        print(f"  factor x{factor} ({label}): " + "  ".join(
            f"{n}={all_metrics[n][f'rmse_ang_{label}']:.2f}" for n in FILTER_SPECS))

    # --- сохранение
    np.savez_compressed(
        f'{RESULTS}/trajectories.npz',
        t=data['t'], q_true=data['q_true'], euler_true=data['euler_true'],
        gyro_bias_true=data['gyro_bias_true'],
        **{f'Q_{n}': trajs[n] for n in trajs},
        **({f'diag_{k}': v for k, v in diag_saved.items()} if diag_saved else {}))
    with open(f'{RESULTS}/metrics_synthetic.json', 'w') as fh:
        json.dump(all_metrics, fh, indent=2)

    # --- markdown-таблицы
    rows_main = ['| Алгоритм | RMSE угла, ° | RMSE roll | RMSE pitch | RMSE yaw '
                 '| Средняя, ° | Макс, ° | мкс/шаг |',
                 '|---|---|---|---|---|---|---|---|']
    for n, m in all_metrics.items():
        rows_main.append(
            f"| {n} | {m['rmse_ang']:.2f} | {m['rmse_roll']:.2f} | "
            f"{m['rmse_pitch']:.2f} | {m['rmse_yaw']:.2f} | {m['mean_ang']:.2f} "
            f"| {m['max_ang']:.2f} | {m['step_us']:.1f} |")

    seg_names = list(SEGMENTS.keys())
    rows_seg = ['| Алгоритм | ' + ' | '.join(seg_names) + ' |',
                '|' + '---|' * (len(seg_names) + 1)]
    for n, m in all_metrics.items():
        rows_seg.append(f"| {n} | " + " | ".join(
            f"{m[f'seg_{s}_ang']:.2f}" for s in seg_names) + " |")

    rows_dt = ['| Алгоритм | 250 Гц | 50 Гц | 21 Гц |', '|---|---|---|---|']
    for n, m in all_metrics.items():
        rows_dt.append(f"| {n} | {m['rmse_ang']:.2f} | "
                       f"{m['rmse_ang_50Hz']:.2f} | {m['rmse_ang_21Hz']:.2f} |")

    with open(f'{RESULTS}/metrics_synthetic.md', 'w', encoding='utf-8') as fh:
        fh.write("# Результаты на синтетическом сценарии (test, seed=42)\n\n"
                 "Все фильтры настроены сеточным поиском на отдельном "
                 "tune-сценарии (seed=1).\nПереходный процесс первых 5 с "
                 "исключён из общих метрик.\n\n## Общие метрики\n\n"
                 + "\n".join(rows_main)
                 + "\n\n## RMSE полного угла по сегментам, °\n\n"
                 + "\n".join(rows_seg)
                 + "\n\n## Устойчивость к частоте дискретизации "
                   "(RMSE угла, °)\n\n" + "\n".join(rows_dt) + "\n\n"
                 "## Параметры после тюнинга\n\n```json\n"
                 + json.dumps({n: best[n]['params'] for n in best}, indent=2)
                 + "\n```\n")
    print(f"\nСохранено: {RESULTS}/metrics_synthetic.md, trajectories.npz")


if __name__ == '__main__':
    main()
