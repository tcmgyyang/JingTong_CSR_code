# -*- coding: utf-8 -*-
"""
Tissue corroboration (GSE153761, cervical cartilage-endplate, 3 degenerate vs 3 normal).
Shows hub-gene direction-of-change in the DISEASE TISSUE, complementing the blood analysis
(blood-silent inflammatory/matrix hub genes are expected to be active here).
"""
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F
F.apply()

ROOT = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.dirname(ROOT)
FIG = os.path.join(ROOT, 'figures'); os.makedirs(FIG, exist_ok=True)
gm = os.path.join(PROJ, '87geoPharm(1)', '04.ann', 'geneMatrix.txt')
expr = pd.read_csv(gm, sep='\t', index_col=0)
groups = {'GSM4653870':'Degenerate','GSM4653871':'Degenerate','GSM4653872':'Degenerate',
          'GSM4653873':'Normal','GSM4653874':'Normal','GSM4653875':'Normal'}
expr = expr[list(groups)]
lab = pd.Series(groups)
deg = [k for k in groups if groups[k]=='Degenerate']; nor=[k for k in groups if groups[k]=='Normal']

hub = ['AKT1','IL1B','IL6','JUN','EGFR','TP53','MMP9','IL4','IL1A','MMP3','VCAM1']
blood_silent = ['JUN','IL4','IL6','VCAM1','IL1A','EGFR','MMP3']  # FPKM<1 in blood (GSE223227)

rows=[]
for h in hub:
    if h not in expr.index:
        rows.append([h, np.nan, np.nan, np.nan, 'not on array']); continue
    v = expr.loc[h]; a=v[deg].astype(float); b=v[nor].astype(float)
    lfc = a.mean()-b.mean()
    try: p = stats.ttest_ind(a,b).pvalue
    except Exception: p = np.nan
    rows.append([h, round(a.mean(),2), round(b.mean(),2), round(lfc,2), round(p,3)])
res = pd.DataFrame(rows, columns=['hub_gene','Degenerate_mean','Normal_mean','log2FC_deg_vs_nor','p']).set_index('hub_gene')
res['blood_expressed'] = ['no (silent)' if h in blood_silent else 'yes' for h in res.index]
res.to_csv(os.path.join(ROOT,'tissue_hubgene_GSE153761.csv'))
print(res.to_string())

# ---- Figure: hub-gene expression in endplate tissue (deg vs normal) ----
# n=3/group: dynamite bar + SD is misleading at this n, so show the three
# individual points per group with a short mean crossbar (ax.hlines) instead.
present = [h for h in hub if h in expr.index]
n = len(present); ncol = 4; nrow = int(np.ceil((n + 1) / ncol))  # +1 cell reserved for the shared legend
ORDER = ['Normal', 'Degenerate']
COL = F.PAL2                                   # {'Normal':'#0072B2','Degenerate':'#D55E00'}
long = expr.loc[present].T.astype(float).copy(); long['group'] = lab
XPOS = {'Normal': 0, 'Degenerate': 1}
RNG = np.random.default_rng(7)                 # deterministic horizontal jitter for the points

fig, axes = plt.subplots(nrow, ncol, figsize=(7.0, 2.35 * nrow))
axes = np.array(axes).ravel()

for k, (ax, h) in enumerate(zip(axes, present)):
    ymax = -np.inf
    for g in ORDER:
        x0 = XPOS[g]; vals = long.loc[lab == g, h].values
        jit = RNG.uniform(-0.13, 0.13, size=vals.size)
        ax.scatter(x0 + jit, vals, s=22, facecolor=COL[g], edgecolor='black',
                   linewidth=0.5, zorder=3, clip_on=False)
        ax.hlines(vals.mean(), x0 - 0.26, x0 + 0.26, color=COL[g], lw=2.0, zorder=2)  # mean crossbar
        ymax = max(ymax, vals.max())
    p = res.loc[h, 'p']; s = F.star(p) if pd.notna(p) else ''
    ax.set_title(f'{h}\n(P = {p}{("  "+s) if s else ""})', fontsize=8)
    ax.set_xlim(-0.55, 1.55); ax.set_xticks([0, 1]); ax.set_xticklabels(ORDER)
    ax.set_xlabel(''); ax.set_ylabel('log2 expression')
    ax.margins(y=0.18)
    F.panel(ax, chr(65 + k))                    # panel letters A, B, C, ...

# blank out any leftover cells, then drop the shared legend into the last cell
for ax in axes[n:]:
    ax.axis('off')
leg_ax = axes[n]                                # first empty cell after the panels
handles = [Line2D([0], [0], marker='o', linestyle='none', markersize=7,
                  markerfacecolor=COL[g], markeredgecolor='black', markeredgewidth=0.5,
                  label=g) for g in ORDER]
handles.append(Line2D([0], [0], color='black', lw=2.0, label='Group mean'))
leg_ax.legend(handles=handles, loc='center', frameon=False, title='Endplate tissue',
              handletextpad=0.6, labelspacing=0.8)
leg_ax.text(0.5, 0.06, 'n = 3 per group', transform=leg_ax.transAxes,
            ha='center', va='center', fontsize=7, color=F.NEUTRAL)

fig.suptitle('Hub-gene expression in cervical endplate tissue (GSE153761)', y=1.005, fontsize=10)
fig.tight_layout()
F.save(fig, os.path.join(FIG, 'Fig_tissue_hubgenes'))

# ---- complementarity summary ----
ts = res.dropna(subset=['log2FC_deg_vs_nor'])
active_silent = ts[(ts.blood_expressed=='no (silent)') & (ts.log2FC_deg_vs_nor.abs()>=0.3)]
print('\nBlood-silent hub genes that ARE differential in endplate tissue (|log2FC|>=0.3):')
print(active_silent[['Degenerate_mean','Normal_mean','log2FC_deg_vs_nor','p']].to_string())
print('\nFigure -> figures/Fig_tissue_hubgenes.png')
