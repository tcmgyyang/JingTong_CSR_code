# -*- coding: utf-8 -*-
"""
Replot of GSVA pathway-activity figure (NO re-run of gseapy gsva).
Reads the saved score matrix + key-pathway stats + sample labels and draws a
publication-ready horizontal grouped box plot.

Inputs (all pre-computed):
    GSVA_kegg_scores.csv   pathway x sample GSVA enrichment-score matrix
    GSVA_key_pathways.csv  pathway, delta_CERVvsCON, p  (network-pharmacology pathways)
    ../ML_GSE223227/labels.csv  sample -> group (CON / CSR / DCM)
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F

HERE = os.path.dirname(os.path.abspath(__file__))
ML   = os.path.join(os.path.dirname(HERE), 'ML_GSE223227')
FIG  = os.path.join(HERE, 'figures')
os.makedirs(FIG, exist_ok=True)
F.apply()

# ---------------- load pre-computed data ----------------
mat    = pd.read_csv(os.path.join(HERE, 'GSVA_kegg_scores.csv'), index_col=0)
key    = pd.read_csv(os.path.join(HERE, 'GSVA_key_pathways.csv'))
labels = pd.read_csv(os.path.join(ML, 'labels.csv'), index_col=0)['group']
mat    = mat[labels.index]                      # column order == label order

# Order pathways by effect size (delta) for a readable, ranked y-axis.
key = key.sort_values('delta_CERVvsCON', ascending=True).reset_index(drop=True)
paths = key['pathway'].tolist()
pmap  = dict(zip(key['pathway'], key['p']))

GROUPS = ['CON', 'CSR', 'DCM']
N = {g: int((labels == g).sum()) for g in GROUPS}
# CB-safe Okabe-Ito groups + a redundant hatch cue so groups separate in grayscale.
HATCH = {'CON': '', 'CSR': '///', 'DCM': '...'}

# ---------------- layout ----------------
# Horizontal: pathways on y-axis so the FULL names fit; 3 boxes per pathway row.
n_path = len(paths)
fig_h  = max(4.2, 0.62 * n_path + 1.4)
fig, ax = plt.subplots(figsize=(7.0, min(fig_h, 7.0)))

row_h   = 1.0                       # vertical slot per pathway
bw      = 0.24                      # box width (in data units along y)
offs    = {'CON': -bw, 'CSR': 0.0, 'DCM': +bw}

for i, path in enumerate(paths):
    y0 = i * row_h
    for g in GROUPS:
        vals = mat.loc[path, labels == g].dropna().values.astype(float)
        if vals.size == 0:
            continue
        bp = ax.boxplot(
            vals, positions=[y0 + offs[g]], widths=bw * 0.92,
            vert=False, patch_artist=True, showfliers=False,
            medianprops=dict(color='black', lw=1.0),
            whiskerprops=dict(color=F.PAL[g], lw=0.9),
            capprops=dict(color=F.PAL[g], lw=0.9),
            boxprops=dict(lw=0.7),
        )
        for box in bp['boxes']:
            box.set_facecolor(F.PAL[g])
            box.set_edgecolor('black')
            box.set_alpha(0.92)
            box.set_hatch(HATCH[g])
        # jittered raw points (small n -> always show the data)
        jit = (np.random.RandomState(i).rand(vals.size) - 0.5) * bw * 0.55
        ax.scatter(vals, np.full(vals.size, y0 + offs[g]) + jit,
                   s=5, color='black', alpha=0.45, linewidths=0, zorder=3)

    # significance star (CSR+DCM vs CON, from the saved p column) at the right margin
    star = F.star(pmap[path])
    xr = mat.loc[path].max()
    t = ax.text(xr + 0.04, y0, star, va='center', ha='left',
                fontsize=9, fontweight='bold', color='black')

# zero reference line (GSVA scores are centred -> 0 = no differential enrichment)
ax.axvline(0, color=F.NEUTRAL, lw=0.9, ls='--', zorder=0)

ax.set_yticks([i * row_h for i in range(n_path)])
ax.set_yticklabels(paths, fontsize=7.5)
ax.set_ylim(-0.7, (n_path - 1) * row_h + 0.7)
ax.set_xlabel('GSVA enrichment score (a.u.)')
ax.set_ylabel('')
ax.margins(x=0.10)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)

# legend: group colours+hatches, placed ABOVE the plot so it never covers data
handles = [mpatches.Patch(facecolor=F.PAL[g], edgecolor='black', hatch=HATCH[g],
                          label=f'{g} (n={N[g]})') for g in GROUPS]
ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.012),
          ncol=3, title=None, fontsize=7.5,
          handlelength=1.6, columnspacing=1.8, borderpad=0.4, frameon=False)
# significance key on its own clean line below the x-axis label
fig.text(0.5, -0.005,
         'CSR+DCM vs CON (Mann-Whitney U): no pathway reached P < 0.05',
         fontsize=6.5, color='#555555', va='top', ha='center')

ax.set_title('GSVA activity of network-pharmacology pathways',
             fontsize=9, pad=34)
fig.tight_layout()
F.save(fig, os.path.join(FIG, 'Fig_GSVA_pathways'))
print('Saved Fig_GSVA_pathways (PNG+PDF). n_path =', n_path, '| groups:', N)
