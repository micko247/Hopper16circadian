#!/usr/bin/env python3
"""heatmap_Cvi_ZT_PIF_SPT.py
-----------------------------
PIF/SPT family — V.vi control, labelled by VIT ID.

VIT_12s0028g01110  Phytochrome-interacting factor 5 PIL6
VIT_07s0005g05100  Phytochrome interacting factor 3-like 5
VIT_07s0005g02510  Phytochrome interacting factor 3 (PIF3)
VIT_14s0060g00260  Phytochrome interacting factor 3 (PIF3)
VIT_18s0001g10270  Basic helix-loop-helix protein SPATULA
VIT_07s0031g00450  Basic helix-loop-helix protein SPATULA

Columns : ZT4, ZT5, ZT6, ZT8, ZT12
          (raw times 24 h, 1 h, 2 h, 4 h, 8 h)

Output  : circadian_heatmap_Cvi_ZT_PIF_SPT.png  (180 dpi)
          circadian_heatmap_Cvi_ZT_PIF_SPT.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ── Expression data — Cabernet Sauvignon control only ────────────────────────
# Internal array order: [1 h, 2 h, 4 h, 8 h, 24 h]
GENES = [
    dict(short='VIT_12s0028g01110', Cab=[13.4242, 13.5538, 13.6000, 12.6323, 14.0519], q_time=9.73e-11, q_culttime=4.20e-8),
    dict(short='VIT_07s0005g05100', Cab=[13.0532, 13.2428, 13.0224, 12.8433, 13.2368], q_time=2.57e-5,  q_culttime=1.99e-3),
    dict(short='VIT_07s0005g02510', Cab=[12.0964, 12.2695, 11.4928, 11.5887, 12.4360], q_time=1.04e-9,  q_culttime=3.83e-4),
    dict(short='VIT_14s0060g00260', Cab=[12.7655, 12.6513, 12.3104, 12.4287, 13.1832], q_time=4.47e-5,  q_culttime=3.23e-8),
    dict(short='VIT_18s0001g10270', Cab=[ 8.3630,  8.7284,  8.8382,  8.7492,  8.9466], q_time=6.68e-4,  q_culttime=4.90e-6),
    dict(short='VIT_07s0031g00450', Cab=[ 6.1315,  6.0866,  6.3291,  6.3190,  5.9455], q_time=3.04e-7,  q_culttime=1.21e-4),
]

# Column reorder: 24 h → ZT4, 1 h → ZT5, 2 h → ZT6, 4 h → ZT8, 8 h → ZT12
ZT_ORDER    = [4, 0, 1, 2, 3]
TIME_LABELS = ['ZT4', 'ZT5', 'ZT6', 'ZT8', 'ZT12']

N_GENES = len(GENES)
N_TIME  = len(TIME_LABELS)

# ── Build z-scored matrix — reordered columns ────────────────────────────────
matrix = np.zeros((N_GENES, N_TIME))
for gi, g in enumerate(GENES):
    vals = np.array(g['Cab'])
    z = (vals - vals.mean()) / vals.std(ddof=1)
    matrix[gi] = z[ZT_ORDER]

# ── Figure geometry ───────────────────────────────────────────────────────────
CELL_H = 0.55
CELL_W = 0.70
L_PAD  = 2.10   # wider to fit VIT ID labels
R_PAD  = 4.20
T_PAD  = 0.90
B_PAD  = 0.45

hm_h = N_GENES * CELL_H
hm_w = N_TIME  * CELL_W
fig_h = hm_h + T_PAD + B_PAD
fig_w = hm_w + L_PAD + R_PAD

fig = plt.figure(figsize=(fig_w, fig_h))

ax = fig.add_axes([L_PAD / fig_w, B_PAD / fig_h, hm_w / fig_w, hm_h / fig_h])

cb_h = hm_h * 0.55
cb_w = 0.18
ax_cb = fig.add_axes([
    (L_PAD + hm_w + 2.90) / fig_w,
    (B_PAD + hm_h * 0.225) / fig_h,
    cb_w / fig_w,
    cb_h / fig_h,
])

# ── Heatmap ───────────────────────────────────────────────────────────────────
vmax = np.abs(matrix).max()
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im = ax.imshow(matrix, aspect='auto', cmap=plt.cm.RdBu_r, norm=norm,
               interpolation='nearest')

# ── X-axis (ZT) — top ────────────────────────────────────────────────────────
ax.set_xticks(range(N_TIME))
ax.set_xticklabels(TIME_LABELS, fontsize=9, fontweight='bold')
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xlabel('Zeitgeber time (control)', fontsize=8.5, labelpad=6)
ax.tick_params(axis='x', length=3)

# ── Y-axis (VIT ID) ───────────────────────────────────────────────────────────
ax.set_yticks(range(N_GENES))
ax.set_yticklabels([g['short'] for g in GENES], fontsize=7.5, family='monospace')
ax.tick_params(axis='y', length=0, pad=4)

# ── Q-value annotations (right) ───────────────────────────────────────────────
def fmt_q(q):
    return f'{q:.1e}' if q < 0.001 else f'{q:.4f}'

ax_right  = (L_PAD + hm_w) / fig_w
ax_bottom = B_PAD / fig_h
ax_height = hm_h  / fig_h

for gi, g in enumerate(GENES):
    y_fig = ax_bottom + ax_height * (1.0 - (gi + 0.5) / N_GENES)
    fig.text(
        ax_right + 0.008, y_fig,
        f"q(Time) = {fmt_q(g['q_time'])}\n"
        f"q(Geno\u00d7Time) = {fmt_q(g['q_culttime'])}",
        va='center', ha='left', fontsize=6.8, color='#333333',
        transform=fig.transFigure,
    )

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar = plt.colorbar(im, cax=ax_cb)
cbar.set_label('z-score', fontsize=8, labelpad=4)
cbar.ax.tick_params(labelsize=7)
cbar.ax.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.6)

# ── Title ─────────────────────────────────────────────────────────────────────
fig.text(
    0.5, 0.995,
    'PIF/SPT family — V.vi control conditions\n'
    'Vitis vinifera cv. Cabernet Sauvignon (Hopper et al. 2016, BMC Plant Biol)\n'
    'Colour: z-score per gene across time-point means  '
    '\u2502  q-values: BH-FDR 3-way ANOVA (full dataset)',
    ha='center', va='top', fontsize=8.2, transform=fig.transFigure,
)

# ── Save ──────────────────────────────────────────────────────────────────────
for fname in ('circadian_heatmap_Cvi_ZT_PIF_SPT.png', 'circadian_heatmap_Cvi_ZT_PIF_SPT.pdf'):
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    print(f'Saved: {fname}')
