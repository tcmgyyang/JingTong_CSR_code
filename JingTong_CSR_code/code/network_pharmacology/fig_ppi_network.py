# -*- coding: utf-8 -*-
"""PPI network of the high-confidence (>=2-source) drug-disease intersection targets,
from STRING (combined score >=0.4). Nodes coloured + sized by degree; the 18 multi-source
hub genes highlighted. Saves Fig_PPI_network + edge/node tables (Cytoscape-ready)."""
import os, sys, urllib.request, urllib.parse
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import networkx as nx
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F; F.apply()
OUT = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(OUT, 'figures')

genes = pd.read_csv(os.path.join(OUT, 'intersection_2source.csv'))['gene'].dropna().astype(str).tolist()
hubs = set(pd.read_csv(os.path.join(OUT, 'hub_genes_2source.csv'))['gene'].astype(str))
print('intersection genes:', len(genes), '| hubs:', len(hubs))

def string_net(gs):
    data = urllib.parse.urlencode({'identifiers': '%0d'.join(gs), 'species': '9606',
                                   'caller_identity': 'jingtong_netpharm'}).encode()
    req = urllib.request.Request('https://string-db.org/api/tsv/network', data=data,
                                 headers={'User-Agent': 'jingtong-netpharm/1.0'})
    txt = urllib.request.urlopen(req, timeout=90).read().decode('utf-8', 'replace')
    e = []
    for ln in txt.splitlines()[1:]:
        f = ln.split('\t')
        if len(f) >= 6: e.append((f[2], f[3], float(f[5])))
    return e

edges = string_net(genes)
G = nx.Graph()
G.add_nodes_from(genes)
for a, b, sc in edges:
    if sc >= 0.4: G.add_edge(a, b, weight=sc)
G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])   # drop isolates
deg = dict(G.degree())
pd.DataFrame([(a, b, round(w['weight'], 3)) for a, b, w in G.edges(data=True)],
             columns=['node1', 'node2', 'combined_score']).to_csv(os.path.join(OUT, 'PPI_edges.csv'), index=False)
pd.DataFrame([(n, deg[n], n in hubs) for n in G.nodes()],
             columns=['gene', 'degree', 'is_hub']).sort_values('degree', ascending=False).to_csv(
    os.path.join(OUT, 'PPI_nodes.csv'), index=False, encoding='utf-8-sig')
print(f'PPI: {G.number_of_nodes()} nodes / {G.number_of_edges()} edges')

# ---- draw: size+colour by degree, hubs outlined/bold ----
fig, ax = plt.subplots(figsize=(7.4, 7.0)); ax.set_aspect('equal')
pos = nx.spring_layout(G, k=1.1/np.sqrt(G.number_of_nodes()), iterations=250, seed=42)
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.18, width=0.4, edge_color='#9AA0A6')
dv = np.array([deg[n] for n in G.nodes()], float)
dn = (dv - dv.min()) / (dv.max() - dv.min() + 1e-9)
sizes = 120 + 1100 * dn
colors = plt.cm.YlOrRd(0.25 + 0.6 * dn)                  # degree gradient (light->deep red-orange)
edgec = ['black' if n in hubs else '#888888' for n in G.nodes()]
lw = [1.3 if n in hubs else 0.4 for n in G.nodes()]
nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=colors,
                       edgecolors=edgec, linewidths=lw)
# labels: hubs bold/larger, others small
nx.draw_networkx_labels(G, pos, ax=ax, labels={n: n for n in G.nodes() if n in hubs},
                        font_size=7.5, font_weight='bold')
nx.draw_networkx_labels(G, pos, ax=ax, labels={n: n for n in G.nodes() if n not in hubs},
                        font_size=5.6, font_color='#333333')
# degree colourbar + hub legend
sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd, norm=plt.Normalize(dv.min(), dv.max()))
cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.01); cb.set_label('Degree (connectivity)', fontsize=7)
cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#FDB863',
                          markeredgecolor='black',markeredgewidth=1.3,markersize=10,label='Hub gene (bold)'),
                   Line2D([0],[0],marker='o',color='none',markerfacecolor='#FDB863',
                          markeredgecolor='#888888',markeredgewidth=0.5,markersize=7,label='Other target')],
          loc='upper left', fontsize=7, frameon=False)
ax.set_title('PPI network of high-confidence drug-disease targets (STRING, score>=0.4)', fontsize=9)
ax.axis('off'); fig.tight_layout()
F.save(fig, os.path.join(FIG, 'Fig_PPI_network'))
print('saved Fig_PPI_network + PPI_edges.csv / PPI_nodes.csv')
