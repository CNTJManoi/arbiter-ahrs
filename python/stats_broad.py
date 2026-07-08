"""
stats_broad.py — статистика по per-trial ошибкам BROAD.

Читает results/broad_report.json (создаётся broad_benchmark.py) и выводит:
парные двусторонние тесты Уилкоксона ARBITER против каждого фильтра,
поправку Холма на пять сравнений, парные ранговые бисериальные корреляции,
оценки Ходжеса–Лемана парных разностей с 95% доверительными интервалами
знаково-ранговой конструкции, а также results/broad_per_trial.csv
(39 испытаний на 6 фильтров).

Запуск: python stats_broad.py
"""
import csv
import json

import numpy as np
from scipy.stats import wilcoxon

FILTERS = ['Complementary', 'Mahony', 'Madgwick', 'MEKF', 'ARBITER', 'VQF']

rep = json.load(open('results/broad_report.json'))
by_trial = {f: rep[f]['tagp']['by_trial'] for f in FILTERS}
trials = sorted(by_trial['ARBITER'])
n = len(trials)


def totals(f):
    return np.array([by_trial[f][t]['total'] for t in trials])


arb = totals('ARBITER')
rows = []
for f in ['VQF', 'Madgwick', 'Complementary', 'Mahony', 'MEKF']:
    d = arb - totals(f)
    _, p = wilcoxon(arb, totals(f))
    ranks = np.argsort(np.argsort(np.abs(d))) + 1.0
    w_plus, w_minus = ranks[d > 0].sum(), ranks[d < 0].sum()
    r_rb = (w_plus - w_minus) / (w_plus + w_minus)
    walsh = np.sort([(d[i] + d[j]) / 2 for i in range(n)
                     for j in range(i, n)])
    hl = float(np.median(walsh))
    mu = n * (n + 1) / 4
    sd = (n * (n + 1) * (2 * n + 1) / 24) ** 0.5
    k = int(np.floor(mu - 1.959964 * sd))
    rows.append([f, p, r_rb, hl, walsh[k], walsh[len(walsh) - 1 - k]])

rows.sort(key=lambda r: r[1])
prev = 0.0
for i, row in enumerate(rows):
    adj = max(prev, min(1.0, (len(rows) - i) * row[1]))
    prev = adj
    row.append(adj)

hdr = f"{'filter':<14}{'p_raw':>11}{'p_holm':>11}{'r_rb':>7}" \
      f"{'HL':>7}{'CI_lo':>8}{'CI_hi':>8}"
print(hdr)
for f, p, r, hl, lo, hi, ph in rows:
    print(f'{f:<14}{p:>11.2e}{ph:>11.2e}{r:>7.2f}{hl:>7.2f}{lo:>8.2f}{hi:>8.2f}')

with open('results/broad_per_trial.csv', 'w', newline='',
          encoding='utf-8') as fh:
    w = csv.writer(fh)
    w.writerow(['trial', 'filter', 'tagp_total_deg', 'tagp_heading_deg',
                'tagp_inclination_deg', 'zero_shot_total_deg'])
    for t in trials:
        for f in FILTERS:
            item = by_trial[f][t]
            zs = rep[f].get('zero_shot', {}).get('by_trial', {}).get(t, {})
            w.writerow([t, f, item.get('total'), item.get('heading'),
                        item.get('inclination'), zs.get('total', '')])
print(f'per-trial CSV: results/broad_per_trial.csv ({n} trials x {len(FILTERS)} filters)')
