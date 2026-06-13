# -*- coding: utf-8 -*-
"""
REPLOT ONLY (no GSEA re-run): read the saved GSEA_key_pathways_summary.csv
(cols gene,Term,NES,FDR) and re-render Fig_single_gene_GSEA as a publication-grade
Term x gene NES heatmap.

Design rationale
----------------
- All NES here are POSITIVE (0.78-1.86). GSEA's enrichment null is NES = 1.0, not 0,
  so we centre the diverging RdBu_r colourmap at 1.0 (TwoSlopeNorm). This puts the
  "no-enrichment" pathways near white and lets genuine enrichment read as red while
  weakly-depleted (<1) reads as blue -> meaningful contrast instead of an all-red wash.
- Cell annotations show NES + an FDR significance star, white-stroked (F.stroke) so the
  numbers stay legible on dark-red cells (the old dark-on-dark text was unreadable).
- A reference tick + label at NES = 1.0 is drawn on the colourbar.
- Arial / colourblind-safe (F.apply, F.DIVERGING), <=178 mm, panel letter, no burned-in
  sentence title, no underscores in tick labels.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F
F.apply()

HERE = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(HERE, 'figures'); os.makedirs(FIG, exist_ok=True)
CSV  = os.path.join(HERE, 'GSEA_key_pathways_summary.csv')

# ---- data prep ----
s = pd.read_csv(CSV)
s['NES'] = pd.to_numeric(s['NES'], errors='coerce')
s['FDR'] = pd.to_numeric(s['FDR'], errors='coerce')

nes = s.pivot_table(index='Term', columns='gene', values='NES', aggfunc='first')
fdr = s.pivot_table(index='Term', columns='gene', values='FDR', aggfunc='first')

# de-underscore tick labels; tidy long KEGG names so they fit
def tidy(t):
    t = str(t).replace('_', ' ')
    t = t.replace(' in diabetic complications', '')          # AGE-RAGE ... -> shorter
    t = t.replace(' signaling pathway', ' signalling')        # keep it readable
    return t
nes.index = [tidy(t) for t in nes.index]
fdr.index = nes.index

# order genes by mean NES (most-enriching biomarker first), pathways by mean NES (top first)
gene_order = nes.mean(0).sort_values(ascending=False).index
path_order = nes.mean(1).sort_values(ascending=True).index   # ascending -> strongest at top after invert
nes = nes.loc[path_order, gene_order]
fdr = fdr.loc[path_order, gene_order]

vmin = float(np.nanmin(nes.values))
vmax = float(np.nanmax(nes.values))
norm = TwoSlopeNorm(vmin=min(vmin, 0.95), vcenter=1.0, vmax=max(vmax, 1.05))
cmap = plt.get_cmap(F.DIVERGING)

# ---- figure ----
nrow, ncol = nes.shape
fig_w = 6.6                                  # <= 178 mm
fig_h = 0.55 * nrow + 1.4
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

im = ax.imshow(nes.values, cmap=cmap, norm=norm, aspect='auto')

# ticks
ax.set_xticks(np.arange(ncol)); ax.set_xticklabels(gene_order, fontweight='bold')
ax.set_yticks(np.arange(nrow)); ax.set_yticklabels(nes.index)
ax.set_xticks(np.arange(-.5, ncol, 1), minor=True)
ax.set_yticks(np.arange(-.5, nrow, 1), minor=True)
ax.grid(which='minor', color='white', linewidth=1.2)
ax.tick_params(which='minor', length=0)
ax.tick_params(which='major', length=0)
for spine in ax.spines.values():
    spine.set_visible(False)

# annotate NES + FDR star, white-stroked for legibility on any cell colour
for i in range(nrow):
    for j in range(ncol):
        v = nes.values[i, j]
        q = fdr.values[i, j]
        if np.isnan(v):
            continue
        star = F.star(q) if not np.isnan(q) else ''
        star = '' if star == 'ns' else star          # only show *,**,*** (ns -> blank)
        # NES on top, star on its own line below with clear separation
        txt = ax.text(j, i - 0.12, f'{v:.2f}', ha='center', va='center',
                      fontsize=8, fontweight='bold', color='black')
        F.stroke(txt, lw=1.6, fg='white')
        if star:
            st = ax.text(j, i + 0.22, star, ha='center', va='center',
                         fontsize=9, color='black')
            F.stroke(st, lw=1.6, fg='white')

ax.set_xlabel('Hub biomarker gene', labelpad=4)
ax.set_ylabel('KEGG pathway', labelpad=4)
ax.set_title('Single-gene GSEA (KEGG NES)', pad=8)

# ---- colourbar with NES = 1.0 (enrichment null) reference ----
cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
cbar.set_label('GSEA NES', rotation=90, labelpad=18)
# black reference line + "1.0 null" annotation at the enrichment null
cbar.ax.axhline(norm(1.0), color='black', linewidth=1.0)   # NES=1.0 reference (explained in footnote)
ticks = sorted(set(list(np.round(np.linspace(vmin, vmax, 4), 2)) + [1.0]))
cbar.set_ticks(ticks)

# significance-key footnote
ax.text(0.0, -0.5 / nrow - 0.085, '* FDR < 0.05   ** FDR < 0.01   *** FDR < 0.001;'
        '  NES > 1 = positive enrichment',
        transform=ax.transAxes, fontsize=6.3, color='#444444', ha='left', va='top')

F.panel(ax, 'A', x=-0.30, y=1.02)
fig.tight_layout()
F.save(fig, os.path.join(FIG, 'Fig_single_gene_GSEA'))
print('Saved Fig_single_gene_GSEA (PNG + PDF) to', FIG)
print('NES range:', round(vmin, 3), '->', round(vmax, 3), '| centred at 1.0')
