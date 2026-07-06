"""
ablation.py — абляционное исследование ARBITER: вклад каждого механизма.

Варианты (каждый выключает ровно один механизм):
  full         — полный алгоритм (TAGP-параметры для BROAD, тюнинг — для синтетики)
  no_budget    — токен-бакеты отключены (cap, rho -> огромные): коррекция
                 всегда полным снапом при открытых воротах
  no_gates     — ворота всегда открыты (пороги гейтов < 0): бюджет остаётся
  no_ub        — самонастройка u_b отключена (rho = rho0 = const)
  no_difftest  — дифференциальные тесты отключены (остаются норма и dip)
  no_dualrate  — burst разрешён всегда (одноуровневый бакет)

Метрики: RMSE полного угла на синтетике (весь сценарий + сегменты помех)
и средний total RMSE по 39 испытаниям BROAD.
Запуск: python ablation.py   (требует results/broad_report.json)
"""
import json
import numpy as np
from multiprocessing import Pool
import os

from simulation import generate_scenario
from benchmark import run_filter, compute_metrics
from broad_benchmark import run_one, field_azimuth, DATA_DIR
from new_filter import ArbiterFilter

DEG = np.pi / 180.0

VARIANTS = {
    'full': {},
    'no_budget': dict(cap_tilt=1e6, cap_yaw=1e6, rho0=1e6),
    'no_gates': dict(gate_hi=-1.0, gate_lo=-2.0),
    'no_ub': dict(ub0=0.0, ub_min=0.0, ub_max=0.0, ub_grow=0.0, ub_boost=0.0),
    'no_difftest': dict(diff_tol=1e6, mag_diff_tol=1e6),
    # отвергнутая ветка: двухскоростная политика расхода (negative result)
    'dualrate': dict(burst_evidence=0.75, routine_mult=2.0),
}


def main():
    with open('results/tuned_params.json') as fh:
        synth_params = json.load(fh)['ARBITER']['params']
    with open('results/broad_report.json', encoding='utf-8') as fh:
        broad_params = json.load(fh)['ARBITER']['tagp']['params']

    data = generate_scenario(seed=42, variant='test')
    trials = sorted(fn for fn in os.listdir(DATA_DIR) if fn.endswith('.hdf5'))
    psi_b = field_azimuth()

    rows = {}
    tasks = []
    for name, over in VARIANTS.items():
        p_syn = dict(synth_params); p_syn.update(over)
        Q, _, _ = run_filter(ArbiterFilter(**p_syn), data)
        m = compute_metrics(Q, data)
        rows[name] = dict(
            synth_all=m['rmse_ang'], synth_mag=m['seg_mag_dist_ang'],
            synth_dyn=m['seg_dyn_acc_ang'], synth_out=m['seg_mag_out_ang'])
        p_br = dict(broad_params); p_br.update(over)
        tasks += [('ARBITER', p_br, t, psi_b, name) for t in trials]

    with Pool(processes=max(2, os.cpu_count() - 2)) as pool:
        res = pool.map(_run, tasks, chunksize=4)
    agg = {}
    for name, total in res:
        agg.setdefault(name, []).append(total)
    for name in rows:
        rows[name]['broad_total'] = float(np.mean(agg[name]))

    lines = ['# Абляционное исследование ARBITER\n',
             '| Вариант | Синтетика: весь сценарий, ° | магн. помеха | '
             'дин. ускорение | без магнитометра | BROAD total (39), ° |',
             '|---|---|---|---|---|---|']
    for name, r in rows.items():
        lines.append(f"| {name} | {r['synth_all']:.2f} | {r['synth_mag']:.2f} "
                     f"| {r['synth_dyn']:.2f} | {r['synth_out']:.2f} | "
                     f"{r['broad_total']:.2f} |")
    with open('results/ablation.md', 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open('results/ablation.json', 'w') as fh:
        json.dump(rows, fh, indent=2)
    print('\n'.join(lines))


def _run(task):
    fname_filter, params, trial, psi_b, variant = task
    _, _, _, m = run_one((fname_filter, params, trial, psi_b))
    return variant, m['total']


if __name__ == '__main__':
    main()
