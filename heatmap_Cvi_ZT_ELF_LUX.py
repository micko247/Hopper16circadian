#!/usr/bin/env python3
"""heatmap_Cvi_ZT_ELF_LUX.py
------------------------------
Evening-complex / photoreceptor genes — both Cabernet Sauvignon and Ramsey control.

Genes (VIT ID, common name):
  VIT_04s0008g00660  ELF3
  VIT_06s0004g06600  ELF4
  VIT_06s0004g05120  LUX
  VIT_05s0077g00940  PHYB
  VIT_02s0025g03530  GSH1

Layout:  two side-by-side heatmap blocks  (Cabernet | Ramsey)
Columns: ZT4, ZT5, ZT6, ZT8, ZT12
         (raw time-points: 24 h, 1 h, 2 h, 4 h, 8 h)

Output:  circadian_heatmap_Cvi_ZT_ELF_LUX.png  (180 dpi)
         circadian_heatmap_Cvi_ZT_ELF_LUX.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm

# ── Data ─────────────────────────────────────────────────────────────────────
# Internal array order from spreadsheet: [1h, 2h, 4h, 8h, 24h]
GENES = [
    dict(vit='VIT_04s0008g00660', name='ELF3',
         Cab=[12.8277, 12.8714, 12.4261, 13.6608, 12.6355],
         Ram=[12.2060, 12.3333, 13.0446, 13.5813, 12.2062],
         q_time=7.88e-18, q_culttime=1.50e-10),
    dict(vit='VIT_06s0004g06600', name='ELF4',
         Cab=[ 5.5653,  5.4036,  5.3573,  6.9768,  5.5619],
         Ram=[ 6.8099,  6.6157,  6.7612,  6.6056,  6.5294],
         q_time=3.86e-7,  q_culttime=4.21e-9),
    dict(vit='VIT_06s0004g05120', name='LUX',
         Cab=[10.6451, 10.4880, 11.5516, 13.9656, 10.1786],
         Ram=[ 9.9162,  9.9633, 12.5525, 13.7847,  9.5528],
         q_time=1.98e-20, q_culttime=3.14e-9),
    dict(vit='VIT_05s0077g00940', name='PHYB',
         Cab=[13.1753, 13.1517, 12.9691, 12.8915, 13.0865],
         Ram=[12.8398, 12.6762, 12.6023, 12.4744, 13.1617],
         q_time=5.48e-4,  q_culttime=2.14e-6),
    dict(vit='VIT_02s0025g03530', name='GSH1',
         Cab=[12.8261, 12.8065, 12.7746, 12.7336, 12.8332],
         Ram=[12.9243, 13.0093, 13.0469, 13.0327, 12.8070],
         q_time=1.43e-5,  q_culttime=3.01e-2),
]

# Column reorder: 24 h → ZT4, 1 h → ZT5, 2 h → ZT6, 4 h → ZT8, 8 h → ZT12
ZT_IDX   = [4, 0, 1, 2, 3]
ZT_LBLS  = ['ZT4', 'ZT5', 'ZT6', 'ZT8', 'ZT12']
N_G = len(GENES)
N_T = len(ZT_LBLS)

# ── Build z-scored matrices ───────────────────────────────────────────────────
def zscore_matrix(key):
    m = np.zeros((N_G, N_T))
    for gi, g in enumerate(GENES):
        v = np.array(g[key])
        z = (v - v.mean()) / v.std(ddof=1)
        m[gi] = z[ZT_IDX]
    return m

mat_cab = zscore_matrix('Cab')
mat_ram = zscore_matrix('Ram')

vmax = max(np.abs(mat_cab).max(), np.abs(mat_ram).max())
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

# ── Figure layout ─────────────────────────────────────────────────────────────
CELL_H = 0.68
CELL_W = 0.76
GAP    = 0.50   # gap between the two blocks
L_PAD  = 2.20   # wide enough for "GENE\nVIT_xxsxxxxgxxxxx"
R_PAD  = 3.30
T_PAD  = 1.10
B_PAD  = 0.45

hm_h = N_G * CELL_H
hm_w = N_T * CELL_W
fig_h = hm_h + T_PAD + B_PAD
fig_w = L_PAD + 2 * hm_w + GAP + R_PAD

fig = plt.figure(figsize=(fig_w, fig_h))

def make_ax(x0_in):
    return fig.add_axes([
        x0_in / fig_w, B_PAD / fig_h,
        hm_w / fig_w,  hm_h / fig_h,
    ])

ax_cab = make_ax(L_PAD)
ax_ram = make_ax(L_PAD + hm_w + GAP)

cb_h  = hm_h * 0.55
cb_w  = 0.18
ax_cb = fig.add_axes([
    (L_PAD + 2*hm_w + GAP + 0.28) / fig_w,
    (B_PAD + hm_h * 0.225) / fig_h,
    cb_w / fig_w,
    cb_h / fig_h,
])

# ── Draw heatmaps ─────────────────────────────────────────────────────────────
def draw_hm(ax, mat, title, show_yticks=True):
    im = ax.imshow(mat, aspect='auto', cmap=plt.cm.RdBu_r, norm=norm,
                   interpolation='nearest')
    ax.set_xticks(range(N_T))
    ax.set_xticklabels(ZT_LBLS, fontsize=9, fontweight='bold')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.tick_params(axis='x', length=3, pad=3)

    if show_yticks:
        ax.set_yticks(range(N_G))
        ax.set_yticklabels(
            [f'{g["name"]}\n{g["vit"]}' for g in GENES],
            fontsize=7.8, linespacing=1.4,
        )
        ax.tick_params(axis='y', length=0, pad=4)
    else:
        ax.set_yticks([])

    ax.set_title(title, fontsize=10, fontweight='bold', pad=10)
    return im

draw_hm(ax_cab, mat_cab, 'Cabernet Sauvignon', show_yticks=True)
im = draw_hm(ax_ram, mat_ram, 'Ramsey',          show_yticks=False)

# ── Q-value annotations (right of Ramsey block) ───────────────────────────────
def fmt_q(q):
    return f'{q:.1e}' if q < 0.001 else f'{q:.4f}'

ram_right  = (L_PAD + 2*hm_w + GAP) / fig_w
ax_bottom  = B_PAD / fig_h
ax_height  = hm_h  / fig_h

for gi, g in enumerate(GENES):
    y_fig = ax_bottom + ax_height * (1.0 - (gi + 0.5) / N_G)
    fig.text(
        ram_right + 0.005, y_fig,
        f"q(Time) = {fmt_q(g['q_time'])}\n"
        f"q(Cult\u00d7Time) = {fmt_q(g['q_culttime'])}",
        va='center', ha='left', fontsize=6.5, color='#333333',
        transform=fig.transFigure,
    )

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar = plt.colorbar(im, cax=ax_cb)
cbar.set_label('z-score', fontsize=8, labelpad=4)
cbar.ax.tick_params(labelsize=7)
cbar.ax.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.6)

# ── Title ─────────────────────────────────────────────────────────────────────
fig.text(
    0.5, 0.998,
    'Evening-complex / photoreceptor genes — control conditions\n'
    'Vitis vinifera (Hopper et al. 2016, BMC Plant Biol)\n'
    'Colour: z-score per gene across time-point means  '
    '\u2502  q-values: BH-FDR 3-way ANOVA (full dataset)',
    ha='center', va='top', fontsize=8.2, transform=fig.transFigure,
)

# ── Save ──────────────────────────────────────────────────────────────────────
for fname in ('circadian_heatmap_Cvi_ZT_ELF_LUX.png',
              'circadian_heatmap_Cvi_ZT_ELF_LUX.pdf'):
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    print(f'Saved: {fname}')
