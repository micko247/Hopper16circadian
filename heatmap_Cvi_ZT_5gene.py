#!/usr/bin/env python3
"""heatmap_Cvi_ZT_5gene.py
--------------------------
As heatmap_Cvi_ZT_4gene.py with GSH1 (VIT_02s0025g03530) added.

Genes   : ELF3a, ELF4a, LUX1a, PHYB, GSH1
Columns : ZT4, ZT5, ZT6, ZT8, ZT12
          (raw times 24 h, 1 h, 2 h, 4 h, 8 h)

Output  : circadian_heatmap_Cvi_ZT_5gene.png  (180 dpi)
          circadian_heatmap_Cvi_ZT_5gene.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ── Expression data — Cabernet Sauvignon control only ────────────────────────
# Internal array order: [1 h, 2 h, 4 h, 8 h, 24 h]
GENES = [
    dict(short='ELF3a', Cab=[12.8277, 12.8714, 12.4261, 13.6608, 12.6355], q_time=7.88e-18, q_culttime=1.50e-10),
    dict(short='ELF4a', Cab=[13.2714, 13.2597, 13.4236, 13.4866, 12.6139], q_time=2.85e-13, q_culttime=7.74e-3),
    dict(short='LUX1a', Cab=[10.6451, 10.4880, 11.5516, 13.9656, 10.1786], q_time=1.98e-20, q_culttime=3.14e-9),
    dict(short='PHYB',  Cab=[13.1753, 13.1517, 12.9691, 12.8915, 13.0865], q_time=5.48e-4,  q_culttime=2.14e-6),
    dict(short='GSH1',  Cab=[12.8261, 12.8065, 12.7746, 12.7336, 12.8332], q_time=1.43e-5,  q_culttime=3.01e-2),
]

# Column reorder: 24 h → ZT4, 1 h → ZT5, 2 h → ZT6, 4 h → ZT8, 8 h → ZT12
ZT_ORDER    = [4, 0, 1, 2, 3]
TIME_LABELS = ['ZT4', 'ZT5', 'ZT6', 'ZT8', 'ZT12']

N_GENES = len(GENES)
N_TIME  = len(TIME_LABELS)

# ── Build z-scored matrix — reordered columns ─────────────────────────────────
matrix = np.zeros((N_GENES, N_TIME))
for gi, g in enumerate(GENES):
    vals = np.array(g['Cab'])
    z = (vals - vals.mean()) / vals.std(ddof=1)
    matrix[gi] = z[ZT_ORDER]

# ── Figure geometry ───────────────────────────────────────────────────────────
CELL_H = 0.55
CELL_W = 0.70
L_PAD  = 1.20
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

# ── Y-axis (gene short name) ──────────────────────────────────────────────────
ax.set_yticks(range(N_GENES))
ax.set_yticklabels([g['short'] for g in GENES], fontsize=9, style='italic')
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
    'Circadian candidate genes — V.vi control conditions\n'
    'Vitis vinifera cv. Cabernet Sauvignon (Hopper et al. 2016, BMC Plant Biol)\n'
    'Colour: z-score per gene across time-point means  '
    '\u2502  q-values: BH-FDR 3-way ANOVA (full dataset)',
    ha='center', va='top', fontsize=8.2, transform=fig.transFigure,
)

# ── Save ──────────────────────────────────────────────────────────────────────
for fname in ('circadian_heatmap_Cvi_ZT_5gene.png', 'circadian_heatmap_Cvi_ZT_5gene.pdf'):
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    print(f'Saved: {fname}')
