"""
broad_benchmark.py — сравнение алгоритмов на реальных данных BROAD
(Laidig, Caruso, Cereatti, Seel 2021, Data 6(7):72, github.com/dlaidig/broad).

Официальный протокол BROAD:
  * метрики: total / heading / inclination RMSE (формулы из example_code/
    broad_utils.py), только по маске movement;
  * TAGP: один общий набор параметров, минимизирующий средний total RMSE
    по всем 39 испытаниям (сеточный поиск) — применяется к каждому фильтру.

Дополнительно наш протокол:
  * zero-shot: параметры, оттюненные на СИНТЕТИЧЕСКОМ сценарии (seed=1),
    без какой-либо подстройки на BROAD — проверка обобщающей способности.

Выравнивание курса: константа-азимут магнитного поля в СК зала (−2.46°,
оценена по ground truth испытания 01) применяется ОДИНАКОВО ко всем
алгоритмам как поворот оценки вокруг вертикали (эквивалент знания
магнитного склонения). Никакого по-испытательного выравнивания нет.

Запуск: python broad_benchmark.py [--quick]   (--quick: только zero-shot)
"""
import json
import os
import sys
import time
from functools import lru_cache
from itertools import product
from multiprocessing import Pool

import numpy as np
import h5py

DATA_DIR = os.path.join('data', 'broad', 'data_hdf5')
RESULTS = 'results'
DEG = np.pi / 180.0

# сокращённые, но охватывающие сетки TAGP (одинаковая процедура для всех)
TAGP_GRIDS = {
    'Complementary': dict(tau_a=[0.5, 1.0, 2.0, 5.0], tau_m=[1.0, 2.0, 5.0, 10.0]),
    'Mahony': dict(kp=[0.5, 1.0, 2.0, 4.0], ki=[0.0, 0.1, 0.3]),
    'Madgwick': dict(beta=[0.01, 0.02, 0.033, 0.05, 0.1, 0.2]),
    'MEKF': dict(sigma_acc_dir=[0.02, 0.05, 0.15], sigma_mag_dir=[0.02, 0.05, 0.2]),
    'ARBITER': dict(rho0=[0.075*DEG, 0.15*DEG, 0.3*DEG, 0.6*DEG],
                    diff_tol=[0.5, 0.8, 1.2],
                    cap_yaw=[5.0*DEG, 15.0*DEG]),
    'VQF': dict(tauAcc=[1.0, 3.0, 9.0], tauMag=[3.0, 9.0, 27.0]),
}

PSI_B = None  # азимут поля, оценивается один раз из trial 01


def make_filter(name, params):
    from filters_baseline import (ComplementaryFilter, MahonyFilter,
                                  MadgwickFilter, MEKF)
    from new_filter import ArbiterFilter
    from vqf_wrapper import VQFWrapper
    cls = dict(Complementary=ComplementaryFilter, Mahony=MahonyFilter,
               Madgwick=MadgwickFilter, MEKF=MEKF, ARBITER=ArbiterFilter,
               VQF=VQFWrapper)[name]
    return cls(**params)


@lru_cache(maxsize=64)
def load_trial(fname):
    path = os.path.join(DATA_DIR, fname)
    with h5py.File(path, 'r') as f:
        acc = f['imu_acc'][:]
        gyr = f['imu_gyr'][:]
        mag = f['imu_mag'][:]
        opt = f['opt_quat'][:]
        mov = f['movement'][:].astype(bool)
        fs = float(f.attrs['sampling_rate'])
    return acc, gyr, mag, opt, mov, fs


def field_azimuth(fname='01_undisturbed_slow_rotation_A.hdf5'):
    """Азимут горизонтальной составляющей поля в СК зала (по ground truth)."""
    acc, gyr, mag, opt, mov, fs = load_trial(fname)
    ok = ~np.isnan(opt[:, 0])
    q = opt[ok]
    mh = mag[ok] / np.linalg.norm(mag[ok], axis=1)[:, None]
    t = 2 * np.cross(q[:, 1:], mh)
    me = mh + q[:, 0][:, None] * t + np.cross(q[:, 1:], t)
    psi = np.arctan2(me[:, 0], me[:, 1])
    return float(np.arctan2(np.mean(np.sin(psi)), np.mean(np.cos(psi))))


def official_metrics(Q_est, opt, mov):
    """Формулы из broad_utils.py: total/heading/inclination RMSE (град)."""
    ok = mov & ~np.isnan(opt[:, 0])
    qe = Q_est[ok] / np.linalg.norm(Q_est[ok], axis=1)[:, None]
    qo = opt[ok] / np.linalg.norm(opt[ok], axis=1)[:, None]
    qo_inv = qo * np.array([1.0, -1.0, -1.0, -1.0])
    w = (qe[:, 0]*qo_inv[:, 0] - qe[:, 1]*qo_inv[:, 1]
         - qe[:, 2]*qo_inv[:, 2] - qe[:, 3]*qo_inv[:, 3])
    z = (qe[:, 0]*qo_inv[:, 3] + qe[:, 1]*qo_inv[:, 2]
         - qe[:, 2]*qo_inv[:, 1] + qe[:, 3]*qo_inv[:, 0])
    total = 2*np.arccos(np.clip(np.abs(w), 0, 1))
    heading = 2*np.arctan(np.abs(z / np.where(np.abs(w) < 1e-12, 1e-12, w)))
    incl = 2*np.arccos(np.clip(np.sqrt(w**2 + z**2), 0, 1))
    return dict(
        total=float(np.degrees(np.sqrt(np.mean(total**2)))),
        heading=float(np.degrees(np.sqrt(np.mean(heading**2)))),
        inclination=float(np.degrees(np.sqrt(np.mean(incl**2)))))


def run_one(task):
    """Задача пула: (filter_name, params_dict, trial_fname) -> метрики."""
    name, params, fname, psi_b = task
    acc, gyr, mag, opt, mov, fs = load_trial(fname)
    dt = 1.0 / fs
    f = make_filter(name, params)
    N = len(acc)
    Q = np.empty((N, 4))
    t0 = time.perf_counter()
    for i in range(N):
        q = f.update((gyr[i, 0], gyr[i, 1], gyr[i, 2]),
                     (acc[i, 0], acc[i, 1], acc[i, 2]),
                     (mag[i, 0], mag[i, 1], mag[i, 2]), dt)
        Q[i, 0], Q[i, 1], Q[i, 2], Q[i, 3] = q
    us = (time.perf_counter() - t0) / N * 1e6
    # константное выравнивание азимута поля (одинаково для всех алгоритмов)
    half = -0.5 * psi_b
    qz = np.array([np.cos(half), 0.0, 0.0, np.sin(half)])
    Qa = np.empty_like(Q)
    Qa[:, 0] = qz[0]*Q[:, 0] - qz[3]*Q[:, 3]
    Qa[:, 1] = qz[0]*Q[:, 1] - qz[3]*Q[:, 2]
    Qa[:, 2] = qz[0]*Q[:, 2] + qz[3]*Q[:, 1]
    Qa[:, 3] = qz[0]*Q[:, 3] + qz[3]*Q[:, 0]
    m = official_metrics(Qa, opt, mov)
    m['us'] = us
    return name, json.dumps(params), fname, m


def main():
    quick = '--quick' in sys.argv
    trials = sorted(fn for fn in os.listdir(DATA_DIR) if fn.endswith('.hdf5'))
    with open(os.path.join(DATA_DIR, 'trials.json')) as fh:
        tinfo = json.load(fh)
    groups = {g['name']: [t + '.hdf5' for t, info in tinfo['trials'].items()
                          if g['name'] in info['groups']]
              for g in tinfo['groups']}
    psi_b = field_azimuth()
    print(f"trials: {len(trials)}, field azimuth psi_B = {psi_b/DEG:.2f} deg")

    with open(f'{RESULTS}/tuned_params.json') as fh:
        synth = json.load(fh)
    # для VQF zero-shot = заводские параметры (его штатный режим out-of-box)
    synth.setdefault('VQF', {'params': {}})

    # -------- режим zero-shot: параметры с синтетики
    tasks = [(n, synth[n]['params'], t, psi_b) for n in TAGP_GRIDS for t in trials]
    # -------- режим TAGP: сетка
    tagp_tasks = []
    if not quick:
        for n, grid in TAGP_GRIDS.items():
            keys = list(grid)
            for combo in product(*(grid[k] for k in keys)):
                p = dict(zip(keys, combo))
                if n == 'ARBITER':
                    p = dict(p)  # cap_tilt по умолчанию (10°)
                tagp_tasks.append((n, p, None, psi_b))
        tagp_tasks = [(n, p, t, psi_b) for (n, p, _, psi) in tagp_tasks
                      for t in trials]
    all_tasks = tasks + tagp_tasks
    print(f"tasks: zero-shot={len(tasks)}, tagp={len(tagp_tasks)}")

    results = {}
    ckpt = f'{RESULTS}/broad_grid_checkpoint.json'
    if os.path.exists(ckpt):
        with open(ckpt) as fh:
            results = json.load(fh)
        done = sum(len(v) for f in results.values() for v in f.values())
        all_tasks = [tk for tk in all_tasks
                     if tk[2] not in results.get(tk[0], {}).get(
                         json.dumps(tk[1]), {})]
        print(f"resume: {done} results from checkpoint, "
              f"{len(all_tasks)} tasks remain")
    t_start = time.time()
    with Pool(processes=max(2, os.cpu_count() - 2)) as pool:
        for k, (name, pjson, fname, m) in enumerate(
                pool.imap_unordered(run_one, all_tasks, chunksize=4)):
            results.setdefault(name, {}).setdefault(pjson, {})[fname] = m
            if (k + 1) % 100 == 0:
                el = time.time() - t_start
                print(f"  {k+1}/{len(all_tasks)} done, {el:.0f}s elapsed",
                      flush=True)
                with open(f'{RESULTS}/broad_grid_checkpoint.json', 'w') as fh:
                    json.dump(results, fh)

    with open(f'{RESULTS}/broad_grid.json', 'w') as fh:
        json.dump(results, fh, indent=1)

    # -------- агрегация
    def mean_metric(per_trial, metric, trial_subset=None):
        names = trial_subset or list(per_trial.keys())
        return float(np.mean([per_trial[t][metric] for t in names]))

    report = {}
    for name in TAGP_GRIDS:
        entry = {}
        # zero-shot
        pj = json.dumps(synth[name]['params'])
        entry['zero_shot'] = dict(
            params=synth[name]['params'],
            total=mean_metric(results[name][pj], 'total'),
            heading=mean_metric(results[name][pj], 'heading'),
            inclination=mean_metric(results[name][pj], 'inclination'),
            us=mean_metric(results[name][pj], 'us'))
        # TAGP
        if not quick:
            best_pj, best_val = None, np.inf
            for pj2, per_trial in results[name].items():
                if len(per_trial) < len(trials):
                    continue
                v = mean_metric(per_trial, 'total')
                if v < best_val:
                    best_val, best_pj = v, pj2
            per_trial = results[name][best_pj]
            entry['tagp'] = dict(
                params=json.loads(best_pj),
                total=best_val,
                heading=mean_metric(per_trial, 'heading'),
                inclination=mean_metric(per_trial, 'inclination'),
                us=mean_metric(per_trial, 'us'),
                by_group={g: dict(
                    total=mean_metric(per_trial, 'total', gt),
                    heading=mean_metric(per_trial, 'heading', gt),
                    inclination=mean_metric(per_trial, 'inclination', gt))
                    for g, gt in groups.items() if gt},
                by_trial=per_trial)
        report[name] = entry

    with open(f'{RESULTS}/broad_report.json', 'w') as fh:
        json.dump(report, fh, indent=2)

    # -------- markdown
    lines = ['# Результаты на реальных данных BROAD (39 испытаний)\n',
             f'Азимут поля в СК зала: {psi_b/DEG:.2f}° (константа, применена '
             f'ко всем алгоритмам одинаково).\n',
             '## TAGP: общий набор параметров, минимизирующий средний total '
             'RMSE по всем испытаниям\n',
             '| Алгоритм | Total RMSE, ° | Heading RMSE, ° | Inclination '
             'RMSE, ° | мкс/шаг | Параметры |', '|---|---|---|---|---|---|']
    if not quick:
        for n in TAGP_GRIDS:
            e = report[n]['tagp']
            lines.append(f"| {n} | {e['total']:.2f} | {e['heading']:.2f} | "
                         f"{e['inclination']:.2f} | {e['us']:.1f} | "
                         f"`{e['params']}` |")
        lines += ['', '### По группам испытаний (total RMSE, °)', '']
        gnames = [g for g in groups if groups[g]]
        lines.append('| Алгоритм | ' + ' | '.join(gnames) + ' |')
        lines.append('|' + '---|'*(len(gnames)+1))
        for n in TAGP_GRIDS:
            bg = report[n]['tagp']['by_group']
            lines.append(f"| {n} | " + " | ".join(
                f"{bg[g]['total']:.2f}" for g in gnames) + " |")
    lines += ['', '## Zero-shot: параметры с синтетического тюнинга '
              '(без подстройки на BROAD)\n',
              '| Алгоритм | Total RMSE, ° | Heading RMSE, ° | Inclination '
              'RMSE, ° |', '|---|---|---|---|']
    for n in TAGP_GRIDS:
        e = report[n]['zero_shot']
        lines.append(f"| {n} | {e['total']:.2f} | {e['heading']:.2f} | "
                     f"{e['inclination']:.2f} |")
    with open(f'{RESULTS}/broad_results.md', 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('Сохранено: results/broad_results.md, broad_report.json')


if __name__ == '__main__':
    main()
