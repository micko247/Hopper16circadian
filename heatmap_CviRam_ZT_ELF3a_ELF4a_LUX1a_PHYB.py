#!/usr/bin/env python3
"""heatmap_CviRam_ZT_ELF3a_ELF4a_LUX1a_PHYB.py
------------------------------------------------
ELF3a, ELF4a, LUX1a, PHYB — PCH1 excluded.

Genotypes: V.vi (top row), Ramsey (bottom row) per gene panel
Columns : ZT4, ZT5, ZT6, ZT8, ZT12
          (raw times 24 h, 1 h, 2 h, 4 h, 8 h)

Output  : circadian_heatmap_CviRam_ZT_ELF3a_ELF4a_LUX1a_PHYB.png  (180 dpi)
          circadian_heatmap_CviRam_ZT_ELF3a_ELF4a_LUX1a_PHYB.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ── Expression data — control only, internal order [1h, 2h, 4h, 8h, 24h] ────
GENES = [
    dict(
        short      = 'ELF3a',
        Cab = [12.8277, 12.8714, 12.4261, 13.6608, 12.6355],
        Ram = [12.2060, 12.3333, 13.0446, 13.5813, 12.2062],
        q_time     = 7.88e-18,
        q_culttime = 1.50e-10,
    ),
    dict(
        short      = 'ELF4a',
        Cab = [13.2714, 13.2597, 13.4236, 13.4866, 12.6139],
        Ram = [12.8437, 12.7323, 13.1760, 13.2619, 12.4391],
        q_time     = 2.85e-13,
        q_culttime = 7.74e-3,
    ),
    dict(
        short      = 'LUX1a',
        Cab = [10.6451, 10.4880, 11.5516, 13.9656, 10.1786],
        Ram = [ 9.9162,  9.9633, 12.5525, 13.7847,  9.5528],
        q_time     = 1.98e-20,
        q_culttime = 3.14e-9,
    ),
    dict(
        short      = 'PHYB',
        Cab = [13.1753, 13.1517, 12.9691, 12.8915, 13.0865],
        Ram = [12.8398, 12.6762, 12.6023, 12.4744, 13.1617],
        q_time     = 5.48e-4,
        q_culttime = 2.14e-6,
    ),
]

GENOTYPES   = ['Cab', 'Ram']
GENO_LABELS = ['V.vi', 'Ramsey']

# Column reorder: 24 h → ZT4, 1 h → ZT5, 2 h → ZT6, 4 h → ZT8, 8 h → ZT12
ZT_ORDER    = [4, 0, 1, 2, 3]
TIME_LABELS = ['ZT4', 'ZT5', 'ZT6', 'ZT8', 'ZT12']

N_GENES = len(GENES)
N_GENO  = len(GENOTYPES)
N_TIME  = len(TIME_LABELS)
N_ROWS  = N_GENES * N_GENO   # 8

# ── Z-scored matrix — per gene across 10 values, columns reordered ────────────
matrix = np.zeros((N_ROWS, N_TIME))
for gi, g in enumerate(GENES):
    all10 = np.concatenate([g[geno] for geno in GENOTYPES])
    mu, sd = all10.mean(), all10.std(ddof=1)
    for ri, geno in enumerate(GENOTYPES):
        z = (np.array(g[geno]) - mu) / sd
        matrix[gi * N_GENO + ri] = z[ZT_ORDER]

# ── Figure geometry ───────────────────────────────────────────────────────────
CELL_H = 0.55
CELL_W = 0.70
L_PAD  = 1.60
R_PAD  = 4.20
T_PAD  = 0.90
B_PAD  = 0.45

hm_h = N_ROWS * CELL_H
hm_w = N_TIME * CELL_W
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

# White separators between gene groups
for gi in range(1, N_GENES):
    ax.axhline(gi * N_GENO - 0.5, color='white', linewidth=2.8, zorder=3)

# ── X-axis (ZT) — top ────────────────────────────────────────────────────────
ax.set_xticks(range(N_TIME))
ax.set_xticklabels(TIME_LABELS, fontsize=9, fontweight='bold')
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xlabel('Zeitgeber time (control)', fontsize=8.5, labelpad=6)
ax.tick_params(axis='x', length=3)

# ── Y-axis (genotype labels) ──────────────────────────────────────────────────
ax.set_yticks(range(N_ROWS))
ax.set_yticklabels(
    [GENO_LABELS[r % N_GENO] for r in range(N_ROWS)],
    fontsize=8, style='italic',
)
ax.tick_params(axis='y', length=0, pad=3)

# ── Gene group labels + q-values ─────────────────────────────────────────────
def fmt_q(q):
    return f'{q:.1e}' if q < 0.001 else f'{q:.4f}'

ax_left   = L_PAD / fig_w
ax_right  = (L_PAD + hm_w) / fig_w
ax_bottom = B_PAD / fig_h
ax_height = hm_h  / fig_h

for gi, g in enumerate(GENES):
    mid_row = gi * N_GENO + (N_GENO - 1) / 2.0
    y_fig = ax_bottom + ax_height * (1.0 - (mid_row + 0.5) / N_ROWS)

    fig.text(
        ax_left - 0.008, y_fig,
        g['short'],
        va='center', ha='right', fontsize=9, style='italic',
        transform=fig.transFigure,
    )

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
    'Circadian candidate genes — V.vi and Ramsey control conditions\n'
    'Vitis spp. leaf transcriptomics (Hopper et al. 2016, BMC Plant Biol)\n'
    'Colour: z-score per gene across genotype\u00d7time means  '
    '\u2502  q-values: BH-FDR 3-way ANOVA (full dataset)',
    ha='center', va='top', fontsize=8.2, transform=fig.transFigure,
)

# ── Save ──────────────────────────────────────────────────────────────────────
for fname in (
    'circadian_heatmap_CviRam_ZT_ELF3a_ELF4a_LUX1a_PHYB.png',
    'circadian_heatmap_CviRam_ZT_ELF3a_ELF4a_LUX1a_PHYB.pdf',
):
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    print(f'Saved: {fname}')
