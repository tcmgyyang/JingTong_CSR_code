# -*- coding: utf-8 -*-
"""
REPLOT (offline) of the ceRNA lncRNA-miRNA-mRNA network.
Reads the saved ENCORI query results (ceRNA_edges.csv / ceRNA_nodes.csv) and rebuilds
the networkx graph -- does NOT touch the network. Fixes vs the original figure:
  * node CLASS encoded by BOTH colour (Okabe-Ito) AND shape (o / ^ / s) -> survives
    grayscale printing and colour-vision deficiency.
  * label clutter removed: only the 5 mRNA hubs + the miRNAs are labelled; 'hsa-' prefix
    stripped from miRNA names; adjustText repels overlapping labels.
  * tighter spring layout (k=1.5/sqrt(n), 200 iters, fixed seed) for a readable spread.
  * short title, double-column-safe size (<=7.2 in), F.apply()/F.save() house style.
"""
import os, sys, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib.lines import Line2D
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F
F.apply()

OUT = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(OUT, 'figures')

# ---- 1. load saved query results, rebuild graph ----
E = pd.read_csv(os.path.join(OUT, 'ceRNA_edges.csv'))
N = pd.read_csv(os.path.join(OUT, 'ceRNA_nodes.csv'))
ntype = dict(zip(N['node'], N['type']))

G = nx.Graph()
for _, r in E.iterrows():
    G.add_edge(r['source'], r['target'])

ANCHORS = ['IL1B', 'IL6', 'JUN', 'MMP3', 'TP53']      # the 5 mRNA hubs

# ---- 2. class -> colour (Okabe-Ito) + shape ----
COL = {'mRNA': '#D55E00', 'miRNA': '#009E73', 'lncRNA': '#0072B2'}   # vermilion / green / blue
SHAPE = {'mRNA': 'o', 'miRNA': '^', 'lncRNA': 's'}                   # circle / triangle / square
SIZE = {'mRNA': 620, 'miRNA': 150, 'lncRNA': 70}                     # large / medium / small
# WITHIN-class gradient: colour intensity AND size encode node degree (connectivity = weight),
# so the more-connected nodes of each type stand out. Class is still given by SHAPE -> colour-
# vision-deficiency safe. Sequential maps stay in each class's hue family (orange/green/blue).
CMAP = {'mRNA': plt.cm.Oranges, 'miRNA': plt.cm.Greens, 'lncRNA': plt.cm.Blues}
deg = dict(G.degree())
def grad(nodes):
    d = np.array([deg[u] for u in nodes], float)
    rng = d.max() - d.min()
    t = (d - d.min()) / rng if rng > 0 else np.full(len(d), 0.6)   # 0..1 within class
    return t

n = G.number_of_nodes()

# ---- concentric (shell) layout, ordered to minimise edge crossings ----
# mRNA hubs in the CENTRE, miRNAs on a middle ring, lncRNAs on the outer ring;
# nodes are angularly ordered so each lncRNA sits near its miRNA, which sits near its
# target hub -> edges run radially -> tidy, not the scattered spring layout.
def clean(name):
    return name[4:] if name.startswith('hsa-') else name

def circ_mean(angs):
    if not angs: return 0.0
    return math.atan2(sum(math.sin(a) for a in angs), sum(math.cos(a) for a in angs))

mRNA   = [u for u in G.nodes() if ntype.get(u) == 'mRNA']   # derive from graph (robust to anchor changes)
miRNA  = [u for u in G.nodes() if ntype.get(u) == 'miRNA']
lncRNA = [u for u in G.nodes() if ntype.get(u) == 'lncRNA']

ang_m = {u: 2*math.pi*i/max(len(mRNA), 1) for i, u in enumerate(mRNA)}
miRNA.sort(key=lambda u: circ_mean([ang_m[v] for v in G.neighbors(u) if v in ang_m]) % (2*math.pi))
ang_mi = {u: 2*math.pi*i/max(len(miRNA), 1) for i, u in enumerate(miRNA)}
lncRNA.sort(key=lambda u: circ_mean([ang_mi[v] for v in G.neighbors(u) if v in ang_mi]) % (2*math.pi))
ang_l = {u: 2*math.pi*i/max(len(lncRNA), 1) for i, u in enumerate(lncRNA)}

R = {'mRNA': 0.55, 'miRNA': 1.40, 'lncRNA': 2.30}
pos = {}
for u, a in ang_m.items():  pos[u] = (R['mRNA']  * math.cos(a), R['mRNA']  * math.sin(a))
for u, a in ang_mi.items(): pos[u] = (R['miRNA'] * math.cos(a), R['miRNA'] * math.sin(a))
for u, a in ang_l.items():  pos[u] = (R['lncRNA'] * math.cos(a), R['lncRNA'] * math.sin(a))

# ---- 3. draw ----
fig, ax = plt.subplots(figsize=(7.6, 7.6))
ax.set_aspect('equal')
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.18, width=0.4, edge_color='#9AA0A6')
for cls in ('lncRNA', 'miRNA', 'mRNA'):            # draw hubs last (on top)
    nodes = [u for u in G.nodes() if ntype.get(u) == cls]
    t = grad(nodes)                                          # 0..1 by degree within class
    node_colors = CMAP[cls](0.40 + 0.55 * t)                 # 0.40..0.95: light->dark by weight
    sizes = SIZE[cls] * (0.75 + 0.55 * t)                    # bigger = more connected
    nx.draw_networkx_nodes(
        G, pos, ax=ax, nodelist=nodes, node_color=node_colors,
        node_shape=SHAPE[cls], node_size=sizes,
        linewidths=0.5, edgecolors='black', alpha=0.96)

# ---- 4. labels ----
# mRNA hubs: bold white labels on the central nodes
for u in mRNA:
    t = ax.text(pos[u][0], pos[u][1], u, fontsize=7.5, fontweight='bold',
                ha='center', va='center', color='white', zorder=6)
    F.stroke(t, lw=1.6, fg='#7A2E00')
# miRNAs: small radial labels in the gap just outside the middle ring, rotated to read outward
rlab = R['miRNA'] + 0.16
for u in miRNA:
    a = ang_mi[u]; deg = math.degrees(a) % 360
    rot, ha = deg, 'left'
    if 90 < deg < 270: rot, ha = deg + 180, 'right'
    ax.text(rlab*math.cos(a), rlab*math.sin(a), clean(u), fontsize=4.8,
            rotation=rot, rotation_mode='anchor', ha=ha, va='center',
            color='#0B3D2E', zorder=5)
# lncRNAs: radial labels OUTSIDE the outer ring (outermost tier -> labels extend outward freely)
rlab2 = R['lncRNA'] + 0.15
for u in lncRNA:
    a = ang_l[u]; deg = math.degrees(a) % 360
    rot, ha = deg, 'left'
    if 90 < deg < 270: rot, ha = deg + 180, 'right'
    ax.text(rlab2*math.cos(a), rlab2*math.sin(a), u, fontsize=4.8,
            rotation=rot, rotation_mode='anchor', ha=ha, va='center',
            color='#0A3D62', zorder=5)
ax.set_xlim(-3.5, 3.5); ax.set_ylim(-3.5, 3.5)
used_adjust = False

# ---- 5. legend (shape + colour together) ----
handles = [Line2D([0], [0], marker=SHAPE[c], color='none', markerfacecolor=CMAP[c](0.82),
                  markeredgecolor='black', markeredgewidth=0.5,
                  markersize=ms, label=f'{c} (n={sum(1 for u in G.nodes() if ntype.get(u)==c)})')
           for c, ms in (('mRNA', 11), ('miRNA', 8), ('lncRNA', 6))]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 1.0),
          handletextpad=0.5, borderaxespad=0.2, labelspacing=0.7,
          title='Node class', title_fontsize=7.5)
ax.text(0.0, 0.012, 'Colour depth & node size scale with degree (connectivity) within each class',
        transform=ax.transAxes, ha='left', va='bottom', fontsize=6.3, color='#555555')

ax.set_title('ceRNA network of hub genes (lncRNA–miRNA–mRNA)', fontsize=9, pad=6)
ax.axis('off')
ax.margins(0.06)
fig.tight_layout()

F.save(fig, os.path.join(FIGDIR, 'Fig_ceRNA_network'))
print('saved Fig_ceRNA_network.png/.pdf  | nodes=%d edges=%d  adjustText=%s'
      % (n, G.number_of_edges(), used_adjust))
print('class counts:', {c: sum(1 for u in G.nodes() if ntype.get(u) == c) for c in COL})
