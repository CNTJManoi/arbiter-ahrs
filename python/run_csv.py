"""
run_csv.py — прогон ARBITER (и базовых фильтров) на реальных данных из CSV.

Формат CSV (заголовок обязателен, разделитель — запятая):
    time,ax,ay,az,gx,gy,gz,mx,my,mz
Единицы:
    time — секунды (монотонно растущее);
    ax..az — м/с^2 (удельная сила: в покое ~(0,0,+9.81) при z вверх);
    gx..gz — рад/с (если у вас °/с — укажите --gyro-deg);
    mx..mz — любые единицы (µT, Гаусс, отсчёты) — используется
             нормированное направление и относительная норма.
    Пропуск магнитометра: пустые поля или NaN.

Примеры:
    python run_csv.py imu_log.csv
    python run_csv.py imu_log.csv --gyro-deg --filter Madgwick --out est.csv

Выход: CSV time,qw,qx,qy,qz,roll,pitch,yaw,Ea,Em,flag_dyn,flag_mag
"""
import argparse
import numpy as np

from imu_math import quat_to_euler_batch
from filters_baseline import (ComplementaryFilter, MahonyFilter,
                              MadgwickFilter, MEKF)
from new_filter import ArbiterFilter

FILTERS = dict(ARBITER=ArbiterFilter, Complementary=ComplementaryFilter,
               Mahony=MahonyFilter, Madgwick=MadgwickFilter, MEKF=MEKF)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--filter', default='ARBITER', choices=list(FILTERS))
    ap.add_argument('--gyro-deg', action='store_true',
                    help='гироскоп в °/с (по умолчанию рад/с)')
    ap.add_argument('--out', default='estimate.csv')
    args = ap.parse_args()

    data = np.genfromtxt(args.csv, delimiter=',', names=True)
    t = data['time']
    gyr = np.stack([data['gx'], data['gy'], data['gz']], axis=1)
    if args.gyro_deg:
        gyr = np.radians(gyr)
    acc = np.stack([data['ax'], data['ay'], data['az']], axis=1)
    mag = np.stack([data['mx'], data['my'], data['mz']], axis=1)
    dt = np.diff(t, prepend=t[0] - (t[1] - t[0]))
    dt = np.clip(dt, 1e-4, 0.5)

    f = FILTERS[args.filter]()
    N = len(t)
    Q = np.empty((N, 4))
    Ea = np.zeros(N); Em = np.zeros(N)
    fd = np.zeros(N); fm = np.zeros(N)
    for i in range(N):
        q = f.update(tuple(gyr[i]), tuple(acc[i]), tuple(mag[i]), float(dt[i]))
        Q[i] = q
        if args.filter == 'ARBITER':
            Ea[i], Em[i] = f.Ea, f.Em
            fd[i], fm[i] = f.flag_dyn, f.flag_mag
    eul = np.degrees(quat_to_euler_batch(Q))
    out = np.column_stack([t, Q, eul, Ea, Em, fd, fm])
    np.savetxt(args.out, out, delimiter=',', fmt='%.6f',
               header='time,qw,qx,qy,qz,roll,pitch,yaw,Ea,Em,flag_dyn,flag_mag',
               comments='')
    print(f'{args.filter}: {N} отсчётов -> {args.out}')


if __name__ == '__main__':
    main()
