# -*- coding: utf-8 -*-
"""PPI topological screening (CytoNCA-style), 3 panels left->right:
   full PPI  --(Degree >= median)-->  sub-network  --(Betweenness >= median)-->  core hubs.
All circular nodes, Reds colormap + size by degree; final hubs in a clean circular layout.
Reuses PPI_edges.csv (>=2-source intersection, STRING >=0.4). Saves Fig_PPI_network."""
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F; F.apply()
OUT = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(OUT, 'figures')

E = pd.read_csv(os.path.join(OUT, 'PPI_edges.csv'))
G = nx.Graph()
for _, r in E.iterrows(): G.add_edge(r['node1'], r['node2'], weight=r['combined_score'])

# ---- stepwise topological filtering ----
def med(d): return float(np.median(list(d.values())))
# both criteria use the FULL-network medians (consistent with the 18-hub definition used downstream)
deg = dict(G.degree()); bet = nx.betweenness_centrality(G)
md, mb = med(deg), med(bet)
keep1 = [n for n in G if deg[n] >= md]                       # panel 2: Degree >= median
G1 = G.subgraph(keep1).copy()
keep2 = [n for n in keep1 if bet[n] >= mb]                   # panel 3: + Betweenness >= median = hubs
G2 = G.subgraph(keep2).copy()
hubs = sorted(G2.nodes(), key=lambda n: -deg[n])
print(f'full {G.number_of_nodes()}n/{G.number_of_edges()}e -> Degree>=median {G1.number_of_nodes()}n '
      f'-> Betweenness>=median {G2.number_of_nodes()}n (hubs)')
print('core hubs:', ', '.join(hubs))
pd.DataFrame([(n, G2.degree(n), round(nx.betweenness_centrality(G2)[n], 4)) for n in hubs],
            columns=['gene','degree','betweenness']).to_csv(os.path.join(OUT,'hub_genes_stepwise.csv'),index=False,encoding='utf-8-sig')

# ---- draw 3 panels (ALL circular/concentric, all nodes are circles) ----
def layout(g, kind):
    if kind == 'circular':
        return nx.circular_layout(g)
    if kind == 'shell':                                   # concentric rings by degree (大圆套小圆)
        order = sorted(g.nodes(), key=lambda n: -g.degree(n))
        n = len(order); s = max(1, round(n / 3))
        shells = [sh for sh in (order[:s], order[s:2*s], order[2*s:]) if sh]   # inner = highest degree
        return nx.shell_layout(g, shells)
    return nx.circular_layout(g)

fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.8))
panels = [(G, 'shell', f'PPI network\n({G.number_of_nodes()} targets, {G.number_of_edges()} edges)', False),
          (G1, 'circular', f'Degree ≥ median\n({G1.number_of_nodes()} nodes)', True),
          (G2, 'circular', f'Betweenness ≥ median\n({G2.number_of_nodes()} core hubs)', True)]
for ax, (g, lay, ttl, lab) in zip(axes, panels):
    pos = layout(g, lay)
    dv = np.array([g.degree(n) for n in g.nodes()], float)
    dn = (dv - dv.min())/(dv.max()-dv.min()+1e-9) if len(dv) and dv.max() > dv.min() else np.full(len(dv), 0.6)
    nx.draw_networkx_edges(g, pos, ax=ax, alpha=0.20, width=0.4, edge_color='#C77')
    nx.draw_networkx_nodes(g, pos, ax=ax, node_shape='o',
                           node_size=(120 + 700*dn) if lab else (40 + 320*dn),
                           node_color=plt.cm.Reds(0.30 + 0.6*dn),
                           edgecolors='#7F1010', linewidths=0.5)
    if lab:
        nx.draw_networkx_labels(g, pos, ax=ax, font_size=6.6, font_weight='bold')
    ax.set_title(ttl, fontsize=8.5); ax.axis('off'); ax.margins(0.12)

# arrows + filter labels between panels (figure coords)
for x, txt in [(0.345, 'Degree ≥ median\n(top 50%)'), (0.675, 'Betweenness ≥ median\n(top 50%)')]:
    fig.patches.append(mpatches.FancyArrow(x-0.03, 0.5, 0.055, 0, transform=fig.transFigure,
                       width=0.006, head_width=0.03, head_length=0.018, color='#555555', length_includes_head=True))
    fig.text(x+0.0, 0.62, txt, ha='center', va='bottom', fontsize=7.5, color='#333333', fontweight='bold')

fig.suptitle('Topological screening of the PPI network for core hub genes', fontsize=10, y=1.02)
fig.tight_layout(rect=[0, 0, 1, 0.98])
F.save(fig, os.path.join(FIG, 'Fig_PPI_network'))
print('saved Fig_PPI_network (3-panel stepwise) + hub_genes_stepwise.csv')
