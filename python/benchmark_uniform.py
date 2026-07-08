"""
benchmark_uniform.py — companion-тест с равномерной дискретизацией.

Тот же тестовый сценарий (seed=42, те же траектория, помехи и модели
сенсоров), но с постоянным шагом 1/250 с и без пропусков данных, чтобы
официальную реализацию VQF можно было включить честно. Пять фильтров
получают те же параметры протокола общего набора (results/tuned_params.json),
VQF работает с авторскими параметрами по умолчанию.

Запуск: python benchmark_uniform.py
Выход:  results/uniform_companion.json + таблица в stdout.
"""
import json

from benchmark import FILTER_SPECS, run_filter, compute_metrics, RESULTS
from simulation import generate_scenario
from vqf_wrapper import VQFWrapper

SEGS = ('dyn_acc', 'vibration', 'mag_dist', 'mag_out', 'combined')


def main():
    with open(f'{RESULTS}/tuned_params.json') as fh:
        best = json.load(fh)
    data = generate_scenario(seed=42, variant='test', uniform=True)
    print(f"uniform companion: N={len(data['dt'])}, dt={data['dt'][0]:.6f} s")

    out = {}
    runners = [(n, cls, best[n]['params']) for n, (cls, _) in
               FILTER_SPECS.items()]
    runners.append(('VQF', VQFWrapper, {}))
    for name, cls, params in runners:
        Q, us, _ = run_filter(cls(**params), data)
        m = compute_metrics(Q, data)
        out[name] = {'rmse_ang': m['rmse_ang'],
                     **{f'seg_{s}_ang': m.get(f'seg_{s}_ang') for s in SEGS}}
        print(f"  {name:14s} total={m['rmse_ang']:6.2f}  " +
              "  ".join(f"{s}={m.get(f'seg_{s}_ang', float('nan')):6.2f}"
                        for s in SEGS))
    with open(f'{RESULTS}/uniform_companion.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('saved -> results/uniform_companion.json')


if __name__ == '__main__':
    main()
