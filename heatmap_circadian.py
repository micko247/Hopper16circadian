#!/usr/bin/env python3
"""heatmap_circadian.py
----------------------
Heatmap of 7 circadian candidate genes in control conditions across two
Vitis genotypes and 5 time points (1, 2, 4, 8, 24 h).

Layout  : rows = genotype (V.vi → V.ri), columns = time (1→24 h)
          7 gene panels separated by white dividers

Data    : Hopper16Supp3.xlsx — sheet 'Hopper.Combined.Leaf.NonParam3w'
          log2-normalised microarray means (n = 3 biological replicates each)

Colour  : z-score per gene, calculated across all 10 values
          (2 genotypes × 5 time points) so temporal and genotypic variation
          are both visible on a common scale

Stats   : BH-FDR q-values from the full 3-way non-parametric ANOVA reported
          in Hopper et al. 2016 BMC Plant Biology 16:118 (GEO: GSE78920)
          q(Time)     = main effect of time across all conditions
          q(Geno×Time)= genotype × time interaction (AdjCultTime)

          NOTE: these q-values were computed on the full dataset (control +
          stress, all genotypes). They are valid for annotating temporal
          significance but do not isolate the control-only time effect.
          For control-specific inference, raw replicates from GSE78920
          should be used.

Output  : circadian_heatmap.png  (180 dpi)
          circadian_heatmap.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ── Expression data extracted from Hopper16Supp3.xlsx ────────────────────────
# Column order within each genotype: 1 h, 2 h, 4 h, 8 h, 24 h (control only)
# Gene order and names as specified; Ramsey excluded.
GENES = [
    dict(
        vit        = 'VIT_04s0008g00660',
        short      = 'ELF3a',
        annotation = 'Early flowering 3',
        Cab = [12.8277, 12.8714, 12.4261, 13.6608, 12.6355],
        Rip = [12.8608, 13.2922, 13.3814, 13.9534, 12.9947],
        q_time     = 7.88e-18,
        q_culttime = 1.50e-10,
    ),
    dict(
        vit        = 'VIT_09s0002g02680',
        short      = 'ELF3b',
        annotation = 'Early flowering 3',
        Cab = [11.5278, 11.3439, 10.7898, 11.8018, 11.1410],
        Rip = [10.9326, 11.4252, 12.0453, 12.7279, 10.9180],
        q_time     = 6.87e-11,
        q_culttime = 5.80e-5,
    ),
    dict(
        vit        = 'VIT_09s0002g01440',
        short      = 'ELF4a',
        annotation = 'Early flowering 4',
        Cab = [13.2714, 13.2597, 13.4236, 13.4866, 12.6139],
        Rip = [13.2908, 13.3402, 13.4791, 13.2043, 12.9921],
        q_time     = 2.85e-13,
        q_culttime = 7.74e-3,
    ),
    dict(
        vit        = 'VIT_06s0004g06600',
        short      = 'LUX4b',
        annotation = 'Unknown protein',
        Cab = [5.5653, 5.4036, 5.3573, 6.9768, 5.5619],
        Rip = [5.6515, 6.2087, 6.3398, 5.9924, 6.2531],
        q_time     = 3.86e-7,
        q_culttime = 4.21e-9,
    ),
    dict(
        vit        = 'VIT_06s0004g05120',
        short      = 'LUX1a',
        annotation = 'ARR1 typeB',
        Cab = [10.6451, 10.4880, 11.5516, 13.9656, 10.1786],
        Rip = [11.6670, 12.5345, 13.2078, 13.9746, 11.3349],
        q_time     = 1.98e-20,
        q_culttime = 3.14e-9,
    ),
    dict(
        vit        = 'VIT_08s0040g00900',
        short      = 'LUX1b',
        annotation = 'Myb family',
        Cab = [12.5113, 12.6271, 12.7406, 12.9790, 12.7318],
        Rip = [13.3819, 13.6313, 13.4512, 13.6088, 13.1151],
        q_time     = 7.85e-9,
        q_culttime = 2.91e-5,
    ),
    dict(
        vit        = 'VIT_05s0077g00940',
        short      = 'PHYB',
        annotation = 'Phytochrome B',
        Cab = [13.1753, 13.1517, 12.9691, 12.8915, 13.0865],
        Rip = [13.3295, 13.2469, 13.1068, 13.0010, 13.0031],
        q_time     = 5.48e-4,
        q_culttime = 2.14e-6,
    ),
]

GENOTYPES   = ['Cab', 'Rip']
GENO_LABELS = ['V.vi', 'V.ri']
TIME_LABELS = ['1 h', '2 h', '4 h', '8 h', '24 h']

N_GENES = len(GENES)
N_GENO  = len(GENOTYPES)
N_TIME  = len(TIME_LABELS)
N_ROWS  = N_GENES * N_GENO   # 14

# ── Build z-scored matrix (14 × 5) ───────────────────────────────────────────
# Each gene is z-scored across its 10 values (2 genotypes × 5 time points).
# This places temporal and genotypic variation on a common scale per gene.
matrix = np.zeros((N_ROWS, N_TIME))

for gi, g in enumerate(GENES):
    all10 = np.concatenate([g[geno] for geno in GENOTYPES])
    mu = all10.mean()
    sd = all10.std(ddof=1)
    for ri, geno in enumerate(GENOTYPES):
        row = gi * N_GENO + ri
        matrix[row] = (np.array(g[geno]) - mu) / sd

# ── Figure geometry (inches) ──────────────────────────────────────────────────
CELL_H = 0.55    # height per heatmap row
CELL_W = 0.70    # width per heatmap column
L_PAD  = 4.20    # left  — gene name + genotype tick labels
R_PAD  = 3.20    # right — q-value annotations
T_PAD  = 0.90    # top   — title + x-axis
B_PAD  = 0.45    # bottom

hm_h = N_ROWS * CELL_H
hm_w = N_TIME * CELL_W
fig_h = hm_h + T_PAD + B_PAD
fig_w = hm_w + L_PAD + R_PAD

fig = plt.figure(figsize=(fig_w, fig_h))

# Heatmap axes
ax = fig.add_axes([
    L_PAD / fig_w,
    B_PAD / fig_h,
    hm_w  / fig_w,
    hm_h  / fig_h,
])

# Colorbar axes (centred vertically beside the heatmap)
cb_h   = hm_h * 0.45
cb_w   = 0.18
cb_gap = 0.20
ax_cb = fig.add_axes([
    (L_PAD + hm_w + cb_gap) / fig_w,
    (B_PAD + hm_h * 0.275)  / fig_h,
    cb_w / fig_w,
    cb_h / fig_h,
])

# ── Draw heatmap ──────────────────────────────────────────────────────────────
vmax = np.abs(matrix).max()
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
cmap = plt.cm.RdBu_r

im = ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm,
               interpolation='nearest')

# ── White separators between gene groups ──────────────────────────────────────
for gi in range(1, N_GENES):
    ax.axhline(gi * N_GENO - 0.5, color='white', linewidth=2.8, zorder=3)

# ── X-axis (time) — top ───────────────────────────────────────────────────────
ax.set_xticks(range(N_TIME))
ax.set_xticklabels(TIME_LABELS, fontsize=9, fontweight='bold')
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xlabel('Time (control, hours after excision)', fontsize=8.5, labelpad=6)
ax.tick_params(axis='x', length=3)

# ── Y-axis (genotype) ─────────────────────────────────────────────────────────
ax.set_yticks(range(N_ROWS))
ax.set_yticklabels(
    [GENO_LABELS[r % N_GENO] for r in range(N_ROWS)],
    fontsize=8,
    style='italic',
)
ax.tick_params(axis='y', length=0, pad=3)

# ── Gene label + q-values (positioned in figure coordinates) ─────────────────
def fmt_q(q):
    """Format q-value: scientific if < 0.001, else 4 d.p."""
    return f'{q:.1e}' if q < 0.001 else f'{q:.4f}'

ax_left   = L_PAD / fig_w
ax_bottom = B_PAD / fig_h
ax_height = hm_h  / fig_h

for gi, g in enumerate(GENES):
    # Fractional row position of the group centre (row 0 = top in imshow)
    mid_row = gi * N_GENO + (N_GENO - 1) / 2.0
    y_fig = ax_bottom + ax_height * (1.0 - (mid_row + 0.5) / N_ROWS)

    # Gene name left of the axes
    fig.text(
        ax_left - 0.008,
        y_fig,
        f"{g['short']}\n{g['vit']}\n{g['annotation']}",
        va='center', ha='right',
        fontsize=7.0,
        fontfamily='monospace',
        transform=fig.transFigure,
    )

    # q-values right of the axes
    ax_right = (L_PAD + hm_w) / fig_w
    fig.text(
        ax_right + 0.008,
        y_fig,
        f"q(Time) = {fmt_q(g['q_time'])}\n"
        f"q(Geno\u00d7Time) = {fmt_q(g['q_culttime'])}",
        va='center', ha='left',
        fontsize=6.8,
        color='#333333',
        transform=fig.transFigure,
    )

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar = plt.colorbar(im, cax=ax_cb)
cbar.set_label('z-score', fontsize=8, labelpad=4)
cbar.ax.tick_params(labelsize=7)
# Dashed zero line
cbar.ax.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.6)

# ── Title ─────────────────────────────────────────────────────────────────────
fig.text(
    0.5, 0.995,
    'Circadian candidate genes — control conditions\n'
    'Vitis spp. leaf transcriptomics (Hopper et al. 2016, BMC Plant Biol)\n'
    'Colour: z-score per gene across all genotype\u00d7time means  '
    '\u2502  q-values: BH-FDR 3-way ANOVA (full dataset)',
    ha='center', va='top',
    fontsize=8.2,
    transform=fig.transFigure,
)

# ── Save ──────────────────────────────────────────────────────────────────────
for fname in ('circadian_heatmap.png', 'circadian_heatmap.pdf'):
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    print(f'Saved: {fname}')
