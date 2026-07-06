"""Объединённая фигура внутренних состояний ARBITER (замена fig4/5/6/7)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from simulation import SEGMENTS

COL_A = '#2a78d6'
COL_M = '#1baf7a'
COL_R = '#e34948'
COL_V = '#4a3aa7'
TRUTH = '#52514e'
SURFACE = '#fcfcfb'
GRID = '#e1e0d9'
MUTED = '#898781'
INK = '#0b0b0b'

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
    'savefig.facecolor': SURFACE, 'axes.edgecolor': '#c3c2b7',
    'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.axisbelow': True, 'font.size': 9.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'figure.dpi': 150,
})

Z = np.load('results/trajectories.npz')
t = Z['t']

fig, axes = plt.subplots(3, 1, figsize=(10, 6.6), sharex=True)

ax = axes[0]
ax.plot(t, Z['diag_Ea'], color=COL_A, lw=1.0, label='E_a (accelerometer)')
ax.plot(t, Z['diag_Em'], color=COL_M, lw=1.0, label='E_m (magnetometer)')
ax.axhline(0.35, color=MUTED, lw=0.8, ls='--')
ax.axhline(0.20, color=MUTED, lw=0.8, ls=':')
ax.set_ylabel('evidence')
ax.legend(ncol=2, fontsize=8.5, loc='lower left')

ax = axes[1]
ax.plot(t, np.degrees(Z['diag_Tt']), color=COL_A, lw=1.0, label='T_t (tilt budget)')
ax.plot(t, np.degrees(Z['diag_Ty']), color=COL_V, lw=1.0, label='T_y (heading budget)')
ax.set_ylabel('budget, deg')
ax.legend(ncol=2, fontsize=8.5, loc='upper left')

ax = axes[2]
ax.plot(t, np.degrees(Z['diag_ub']), color=COL_R, lw=1.0,
        label='u_b (drift uncertainty)')
ax.fill_between(t, 0, Z['diag_static'] * np.degrees(Z['diag_ub']).max(),
                color=COL_A, alpha=0.18, lw=0, label='rest detected')
ax.set_ylabel('u_b, deg/s')
ax.set_xlabel('time, s')
ax.legend(ncol=2, fontsize=8.5, loc='upper left')

marks = {'dyn_acc': 'accel.\ntransients', 'vibration': 'vibration',
         'mag_dist': 'magnetic\ndisturbance', 'mag_out': 'magnetometer\noutage'}
for k, ax in enumerate(axes):
    for name, (t0, t1) in SEGMENTS.items():
        if name in marks:
            ax.axvspan(t0, t1, color='#0b0b0b', alpha=0.05, lw=0)
            if k == 0:
                ax.text((t0 + t1) / 2, 1.04, marks[name], ha='center',
                        va='bottom', fontsize=7.5, color=MUTED)

fig.align_ylabels(axes)
fig.tight_layout()
fig.savefig('plots/fig11_internal.png', bbox_inches='tight')
print('plots/fig11_internal.png')
