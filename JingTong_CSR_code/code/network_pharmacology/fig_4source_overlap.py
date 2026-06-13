# -*- coding: utf-8 -*-
"""Four-source drug-target overlap figures: (1) 4-set Venn, (2) UpSet plot.
Shows how TCMSP / STITCH / SwissTargetPrediction / HERB contribute & overlap for the
7 Jingtong-Granules herbs. Reads drug_targets_4source.csv (source-annotated)."""
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F
F.apply()
OUT = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(OUT, 'figures')

d = pd.read_csv(os.path.join(OUT, 'drug_targets_4source.csv'))
SRC = ['TCMSP', 'STITCH', 'SwissTarget', 'HERB']
sets = {s: set(d.loc[d[s] == True, 'gene']) for s in SRC}
print('set sizes:', {s: len(v) for s, v in sets.items()})

# ---- (1) 4-set Venn (Okabe-Ito) ----
try:
    import venn
    COL = ['#0072B2', '#D55E00', '#009E73', '#E69F00']     # Okabe-Ito for the 4 sources
    ax = venn.venn(sets, cmap=COL, fontsize=8, legend_loc='upper left')
    fig = ax.figure; fig.set_size_inches(7.2, 6.2)
    fig.suptitle('Jingtong-Granules drug targets across 4 databases (7 herbs)', fontsize=10, y=0.99)
    F.save(fig, os.path.join(FIG, 'Fig_4source_venn'))
    print('saved Fig_4source_venn')
except Exception as e:
    print('venn failed:', e)

# ---- (2) UpSet plot (clearest for >=4 sets) ----
try:
    from upsetplot import UpSet, from_contents
    data = from_contents(sets)
    fig = plt.figure(figsize=(7.4, 4.6))
    up = UpSet(data, sort_by='cardinality', sort_categories_by='cardinality',
               facecolor='#0072B2', shading_color='#EFEFEF', element_size=None,
               show_counts=True, min_subset_size=2)
    up.plot(fig=fig)
    fig.suptitle('Drug-target overlap across TCMSP / STITCH / SwissTarget / HERB', fontsize=10)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, 'Fig_4source_upset.' + ext), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved Fig_4source_upset')
except Exception as e:
    print('upset failed:', e)

# stats
print('\nunion:', len(set().union(*sets.values())),
      '| >=2 sources:', int((d.n_sources >= 2).sum()),
      '| all-4-shared:', int((d.n_sources == 4).sum()))
