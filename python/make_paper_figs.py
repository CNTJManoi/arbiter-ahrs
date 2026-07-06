"""Генерация языко-зависимых комплектов рисунков статьи.

Выход: plots/en/fig_{pipeline,error,internal,rmse,broad}.png
       plots/ru/fig_{pipeline,error,internal,rmse,broad}.png
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

from simulation import SEGMENTS
from imu_math import quat_angle_batch

COL = {'Complementary': '#2a78d6', 'Mahony': '#1baf7a', 'Madgwick': '#eda100',
       'MEKF': '#4a3aa7', 'ARBITER': '#e34948', 'VQF': '#008300'}
TRUTH = '#52514e'
SURFACE = '#fcfcfb'
GRID = '#e1e0d9'
MUTED = '#898781'
INK = '#0b0b0b'
FILTERS = ['Complementary', 'Mahony', 'Madgwick', 'MEKF', 'ARBITER']

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
    'savefig.facecolor': SURFACE, 'axes.edgecolor': '#c3c2b7',
    'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.axisbelow': True, 'font.size': 9.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'figure.dpi': 150,
    # требования IEEE к графике: шрифт из рекомендованного списка,
    # math в стиле Times; PDF с внедрёнными (Type 42) шрифтами
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'stix',
    'pdf.fonttype': 42,
})


def save_fig(fig, out_png):
    """PNG 300 dpi (для DOCX) + векторный PDF (для подачи в IEEE)."""
    fig.savefig(out_png, bbox_inches='tight', dpi=300)
    fig.savefig(out_png[:-4] + '.pdf', bbox_inches='tight')
    plt.close(fig)

DISPLAY = {
    'en': {},
    'ru': {'Complementary': 'Комплементарный'},
}

STR = {
    'en': dict(
        seg=dict(dyn_acc='accel.\ntransients', vibration='vibration',
                 mag_dist='magnetic\ndisturbance', mag_out='magnetometer\noutage'),
        err_y='orientation error, deg', time='time, s',
        ev='evidence', budget='budget, deg', ub_y='$u_b$, deg/s',
        ev_a='$E_a$ (accelerometer)', ev_m='$E_m$ (magnetometer)',
        b_t='$T_t$ (tilt budget)', b_y='$T_y$ (heading budget)',
        ub='$u_b$ (drift uncertainty)', rest='rest detected',
        rmse_y='total-angle RMSE, deg',
        cases=['whole\nscenario', 'accel.\ntransients', 'vibration',
               'magnetic\ndisturbance', 'magnetometer\noutage', 'combined'],
        broad_cats=['Total RMSE', 'Inclination RMSE', 'Heading RMSE'],
        broad_y='RMSE, deg',
        pipe=dict(
            gyro='Gyroscope\n' + r'$g = \omega + b$' + '\n(rad/s)',
            acc='Accelerometer\n' + r'$a$ (m/s$^2$)',
            mag='Magnetometer\n' + r'$m$',
            pred='Kinematic\nprediction\n' + r'$q \leftarrow q \otimes \exp(\frac{1}{2}\hat\omega \Delta t)$',
            tests_a='Accelerometer tests:\nnorm; differential\n(transport residual)',
            tests_m='Magnetometer tests:\nnorm; dip; differential\n(references learned online)',
            ev_a='Evidence $E_a$\nslow charge, fast discharge\n+ hysteresis gate',
            ev_m='Evidence $E_m$\n(slow recovery)\n+ hysteresis gate',
            bucket='Token buckets\n' + r'$T \leftarrow \min(T + \rho\,\Delta t,\ C)$' + '\n'
                   + r'$\rho = \rho_0 + u_b$',
            tilt='Tilt correction\ngeodesic step\n' + r'$\theta_c = \min(\theta, T_t)$',
            yaw='Heading correction\nabout the vertical\n' + r'$\psi_c = \pm\min(|\psi|, T_y)$',
            rest='Rest detector\nbias by averaging,\nuncertainty decays',
            out='$\\hat q$,\nRPY',
            sat='saturation\n' + r'$\to u_b\uparrow$',
            lbl_b=r'$\hat b$', lbl_ub_dn=r'$u_b\downarrow$',
            lbl_Tt=r'$T_t$', lbl_Ty=r'$T_y$')),
    'ru': dict(
        seg=dict(dyn_acc='броски\nускорения', vibration='вибрация',
                 mag_dist='магнитная\nпомеха', mag_out='магнитометр\nотключён'),
        err_y='ошибка ориентации, °', time='время, с',
        ev='доказательства', budget='бюджет, °', ub_y='$u_b$, °/с',
        ev_a='$E_a$ (акселерометр)', ev_m='$E_m$ (магнитометр)',
        b_t='$T_t$ (бюджет наклона)', b_y='$T_y$ (бюджет курса)',
        ub='$u_b$ (неопределённость дрейфа)', rest='обнаружен покой',
        rmse_y='RMSE полного угла, °',
        cases=['весь\nсценарий', 'броски\nускорения', 'вибрация',
               'магнитная\nпомеха', 'без\nмагнитометра', 'комбин.'],
        broad_cats=['Полный RMSE', 'RMSE наклона', 'RMSE курса'],
        broad_y='RMSE, °',
        pipe=dict(
            gyro='Гироскоп\n' + r'$g = \omega + b$' + '\n(рад/с)',
            acc='Акселерометр\n' + r'$a$ (м/с$^2$)',
            mag='Магнитометр\n' + r'$m$',
            pred='Кинематическое\nпредсказание\n' + r'$q \leftarrow q \otimes \exp(\frac{1}{2}\hat\omega \Delta t)$',
            tests_a='Тесты акселерометра:\nнорма; дифференциальный\n(невязка переноса)',
            tests_m='Тесты магнитометра:\nнорма; наклонение;\nдифференциальный\n(эталоны обучаются)',
            ev_a='Накопитель $E_a$\nмедленный заряд,\nбыстрый разряд\n+ гистерезисные ворота',
            ev_m='Накопитель $E_m$\n(медленный возврат доверия)\n+ гистерезисные ворота',
            bucket='Токен-бакеты\n' + r'$T \leftarrow \min(T + \rho\,\Delta t,\ C)$' + '\n'
                   + r'$\rho = \rho_0 + u_b$',
            tilt='Коррекция наклона\nгеодезический шаг\n' + r'$\theta_c = \min(\theta, T_t)$',
            yaw='Коррекция курса\nвокруг вертикали\n' + r'$\psi_c = \pm\min(|\psi|, T_y)$',
            rest='Детектор покоя\nсмещение — усреднением,\nнеопределённость спадает',
            out='$\\hat q$,\nRPY',
            sat='насыщение\n' + r'$\to u_b\uparrow$',
            lbl_b=r'$\hat b$', lbl_ub_dn=r'$u_b\downarrow$',
            lbl_Tt=r'$T_t$', lbl_Ty=r'$T_y$')),
}

Z = np.load('results/trajectories.npz')
t = Z['t']
q_true = Z['q_true']
with open('results/broad_report.json', encoding='utf-8') as fh:
    BROAD = json.load(fh)


def shade(ax, L, labels=False, ytop=None):
    for name, (t0, t1) in SEGMENTS.items():
        if name in L['seg']:
            ax.axvspan(t0, t1, color='#0b0b0b', alpha=0.05, lw=0)
            if labels:
                yy = ytop if ytop is not None else ax.get_ylim()[1]
                ax.text((t0+t1)/2, yy, L['seg'][name], ha='center', va='top',
                        fontsize=7.5, color=MUTED)


def fig_error(L, out, disp):
    fig, ax = plt.subplots(figsize=(10, 3.6))
    for n in FILTERS:
        err = np.degrees(quat_angle_batch(Z[f'Q_{n}'], q_true))
        ax.semilogy(t, np.maximum(err, 1e-3), color=COL[n], lw=1.1,
                    label=disp.get(n, n))
    ax.set_xlabel(L['time'])
    ax.set_ylabel(L['err_y'])
    ax.set_ylim(0.02, 300)
    shade(ax, L, labels=True, ytop=250)
    ax.legend(ncol=5, loc='lower left', bbox_to_anchor=(0, 1.02), fontsize=8.5)
    fig.tight_layout()
    save_fig(fig, out)


def fig_internal(L, out):
    fig, axes = plt.subplots(3, 1, figsize=(10, 6.6), sharex=True)
    ax = axes[0]
    ax.plot(t, Z['diag_Ea'], color=COL['Complementary'], lw=1.0, label=L['ev_a'])
    ax.plot(t, Z['diag_Em'], color=COL['Mahony'], lw=1.0, label=L['ev_m'])
    ax.axhline(0.35, color=MUTED, lw=0.8, ls='--')
    ax.axhline(0.20, color=MUTED, lw=0.8, ls=':')
    ax.set_ylabel(L['ev'])
    ax.legend(ncol=2, fontsize=8.5, loc='lower left')
    ax = axes[1]
    ax.plot(t, np.degrees(Z['diag_Tt']), color=COL['Complementary'], lw=1.0,
            label=L['b_t'])
    ax.plot(t, np.degrees(Z['diag_Ty']), color=COL['MEKF'], lw=1.0,
            label=L['b_y'])
    ax.set_ylabel(L['budget'])
    ax.legend(ncol=2, fontsize=8.5, loc='upper left')
    ax = axes[2]
    ax.plot(t, np.degrees(Z['diag_ub']), color=COL['ARBITER'], lw=1.0,
            label=L['ub'])
    ax.fill_between(t, 0, Z['diag_static'] * np.degrees(Z['diag_ub']).max(),
                    color=COL['Complementary'], alpha=0.18, lw=0,
                    label=L['rest'])
    ax.set_ylabel(L['ub_y'])
    ax.set_xlabel(L['time'])
    ax.legend(ncol=2, fontsize=8.5, loc='upper left')
    for k, ax in enumerate(axes):
        for name, (t0, t1) in SEGMENTS.items():
            if name in L['seg']:
                ax.axvspan(t0, t1, color='#0b0b0b', alpha=0.05, lw=0)
                if k == 0:
                    ax.text((t0+t1)/2, 1.04, L['seg'][name], ha='center',
                            va='bottom', fontsize=7.5, color=MUTED)
    fig.align_ylabels(axes)
    fig.tight_layout()
    save_fig(fig, out)


def fig_rmse(L, out, disp):
    with open('results/metrics_synthetic.json', encoding='utf-8') as fh:
        M = json.load(fh)
    keys = ['rmse_ang', 'seg_dyn_acc_ang', 'seg_vibration_ang',
            'seg_mag_dist_ang', 'seg_mag_out_ang', 'seg_combined_ang']
    fig, ax = plt.subplots(figsize=(9.5, 3.6))
    x = np.arange(len(keys))
    w = 0.15
    for i, n in enumerate(FILTERS):
        vals = [M[n][k] for k in keys]
        ax.bar(x + (i-2)*w, vals, w*0.92, color=COL[n], label=disp.get(n, n))
    ax.set_yscale('log')
    ax.set_xticks(x, L['cases'], fontsize=8.5)
    ax.set_ylabel(L['rmse_y'])
    ax.legend(ncol=5, fontsize=8.5, loc='lower left', bbox_to_anchor=(0, 1.0))
    ax.grid(axis='x', visible=False)
    fig.tight_layout()
    save_fig(fig, out)


def fig_broad(L, out, disp):
    bf = [n for n in COL if n in BROAD and 'tagp' in BROAD[n]]
    fig, ax = plt.subplots(figsize=(9.5, 3.4))
    cats = ['total', 'inclination', 'heading']
    x = np.arange(len(cats))
    w6 = 0.8 / len(bf)
    for i, n in enumerate(bf):
        vals = [BROAD[n]['tagp'][c] for c in cats]
        bars = ax.bar(x + (i - len(bf)/2 + 0.5)*w6, vals, w6*0.9,
                      color=COL[n], label=disp.get(n, n))
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width()/2, v + 0.08, f'{v:.1f}',
                    ha='center', fontsize=7, color=INK)
    ax.set_xticks(x, L['broad_cats'])
    ax.set_ylabel(L['broad_y'])
    ax.legend(ncol=6, fontsize=8.5, loc='lower left', bbox_to_anchor=(0, 1.0))
    ax.grid(axis='x', visible=False)
    fig.tight_layout()
    save_fig(fig, out)


def fig_pipeline(L, out):
    """Схема: «хребет» кватерниона слева направо, сенсорные ветки строго под
    своими коррекциями, бюджеты сверху, покой слева сверху, петля насыщения —
    пунктиром. Все стрелки ортогональны, пересечений нет."""
    P = L['pipe']
    fig, ax = plt.subplots(figsize=(11.5, 6.0))
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    def box(x0, y0, x1, y1, text, fc='#ffffff', fs=8.3):
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                    boxstyle='round,pad=0.008',
                                    fc=fc, ec='#b7b5ac', lw=1.1))
        ax.text((x0 + x1) / 2, (y0 + y1) / 2, text, ha='center', va='center',
                fontsize=fs)

    def arrow(x1, y1, x2, y2, dashed=False):
        ax.add_patch(FancyArrowPatch(
            (x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=10,
            color=MUTED, lw=1.2, linestyle=(0, (4, 3)) if dashed else 'solid',
            shrinkA=0, shrinkB=0))

    def elbow(pts, dashed=False):
        for a, b in zip(pts[:-2], pts[1:-1]):
            ax.plot([a[0], b[0]], [a[1], b[1]], color=MUTED, lw=1.2,
                    linestyle=(0, (4, 3)) if dashed else 'solid')
        arrow(*pts[-2], *pts[-1], dashed=dashed)

    def lab(x, y, text, ha='left', va='center', fs=7.6):
        ax.text(x, y, text, fontsize=fs, color=INK, ha=ha, va=va)

    # ---- пояса: top 0.78..0.99, spine 0.52..0.72, acc 0.27..0.46, mag 0.02..0.21
    box(0.02, 0.52, 0.14, 0.72, P['gyro'], fc='#eef4fc')
    box(0.20, 0.52, 0.38, 0.72, P['pred'])
    box(0.50, 0.52, 0.68, 0.72, P['tilt'], fc='#f2effc')
    box(0.72, 0.52, 0.90, 0.72, P['yaw'], fc='#f2effc')
    box(0.925, 0.53, 0.99, 0.71, P['out'], fc='#eef4fc')
    box(0.20, 0.78, 0.38, 0.99, P['rest'], fc='#eef9f3')
    box(0.50, 0.78, 0.90, 0.99, P['bucket'], fc='#eef9f3')
    box(0.02, 0.27, 0.14, 0.46, P['acc'], fc='#eef4fc')
    box(0.20, 0.27, 0.44, 0.46, P['tests_a'], fc='#fdf6e8')
    box(0.50, 0.27, 0.68, 0.46, P['ev_a'], fc='#fdecec')
    box(0.02, 0.02, 0.14, 0.21, P['mag'], fc='#eef4fc')
    box(0.20, 0.02, 0.44, 0.21, P['tests_m'], fc='#fdf6e8')
    box(0.72, 0.02, 0.90, 0.21, P['ev_m'], fc='#fdecec')

    # ---- хребет
    arrow(0.14, 0.62, 0.20, 0.62)
    arrow(0.38, 0.62, 0.50, 0.62)
    arrow(0.68, 0.62, 0.72, 0.62)
    arrow(0.90, 0.62, 0.925, 0.62)
    # ---- ветка акселерометра
    arrow(0.14, 0.365, 0.20, 0.365)
    arrow(0.44, 0.365, 0.50, 0.365)
    arrow(0.59, 0.46, 0.59, 0.52)
    # ---- ветка магнитометра
    arrow(0.14, 0.115, 0.20, 0.115)
    arrow(0.44, 0.115, 0.72, 0.115)
    arrow(0.81, 0.21, 0.81, 0.52)
    # ---- покой: вход от гироскопа, выходы b̂ и u_b
    elbow([(0.08, 0.72), (0.08, 0.885), (0.20, 0.885)])
    arrow(0.29, 0.78, 0.29, 0.72)
    lab(0.30, 0.755, P['lbl_b'])
    arrow(0.38, 0.925, 0.50, 0.925)
    lab(0.44, 0.945, P['lbl_ub_dn'], ha='center')
    # ---- бюджеты в коррекции
    arrow(0.55, 0.78, 0.55, 0.72)
    lab(0.535, 0.75, P['lbl_Tt'], ha='right')
    arrow(0.77, 0.78, 0.77, 0.72)
    lab(0.755, 0.75, P['lbl_Ty'], ha='right')
    # ---- петля насыщения (пунктир вверх)
    arrow(0.63, 0.72, 0.63, 0.78, dashed=True)
    arrow(0.85, 0.72, 0.85, 0.78, dashed=True)
    lab(0.865, 0.75, P['sat'], fs=7.4)
    save_fig(fig, out)


if __name__ == '__main__':
    for lang in ('en', 'ru'):
        L = STR[lang]
        disp = DISPLAY[lang]
        d = f'plots/{lang}'
        os.makedirs(d, exist_ok=True)
        fig_pipeline(L, f'{d}/fig_pipeline.png')
        fig_error(L, f'{d}/fig_error.png', disp)
        fig_internal(L, f'{d}/fig_internal.png')
        fig_rmse(L, f'{d}/fig_rmse.png', disp)
        fig_broad(L, f'{d}/fig_broad.png', disp)
        print(lang, 'figures done')
