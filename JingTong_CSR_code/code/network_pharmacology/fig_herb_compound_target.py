# -*- coding: utf-8 -*-
"""Herb-compound-target interaction network: 7 Jingtong-Granules herbs -> active
compounds -> shared (drug-disease intersection) targets. Built from the compound-resolved
sources (TCMSP + SwissTarget + STITCH), restricted to the >=2-source intersection targets."""
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F; F.apply()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(OUT, 'figures')
U = lambda s: str(s).strip().upper()
HERB_CN = {'sanqi': 'Sanqi', 'chuanxiong': 'Chuanxiong', 'yanhusuo': 'Yanhusuo',
           'qianghuo': 'Qianghuo', 'baishao': 'Baishao', 'weilingxian': 'Weilingxian', 'gegen': 'Gegen'}

targets = set(pd.read_csv(os.path.join(OUT, 'intersection_2source.csv'))['gene'].dropna().map(U))
hubs = set(pd.read_csv(os.path.join(OUT, 'hub_genes_2source.csv'))['gene'].map(U))

hc = set()   # (herb, compound)
ct = set()   # (compound, target)
# TCMSP
t = pd.read_csv(os.path.join(ROOT, 'allTargets.symbol.xls'), sep='\t')
t = t[t.Drug.isin(HERB_CN)]
for _, r in t.iterrows():
    g = U(r['Symbol'])
    if g in targets:
        c = str(r['MolName']); hc.add((r['Drug'], c)); ct.add((c, g))
# SwissTarget (herb, gene, compound)
p = os.path.join(OUT, 'swiss_targets.csv')
if os.path.exists(p):
    for _, r in pd.read_csv(p).iterrows():
        g = U(r['gene'])
        if g in targets:
            c = str(r['compound']); hc.add((r['herb'], c)); ct.add((c, g))
# STITCH (herb, gene, MolName)
p = os.path.join(OUT, 'STITCH_targets.csv')
if os.path.exists(p):
    for _, r in pd.read_csv(p).iterrows():
        g = U(r['gene'])
        if g in targets and pd.notna(r.get('MolName')):
            c = str(r['MolName']); hc.add((r['herb'], c)); ct.add((c, g))

compounds = sorted({c for _, c in hc} | {c for c, _ in ct})
print(f'herbs 7 | compounds {len(compounds)} | targets {len(targets)} | hc {len(hc)} ct {len(ct)}')

G = nx.Graph()
for h, c in hc: G.add_edge('H:' + h, 'C:' + c)
for c, g in ct: G.add_edge('C:' + c, 'T:' + g)
# keep only compounds connected to a kept target (drop dangling)
G.remove_nodes_from([n for n in list(G) if n.startswith('C:') and not any(nb.startswith('T:') for nb in G.neighbors(n))])
deg = dict(G.degree())

import math
def typ(n): return n[0]
SHAPE = {'H': 'H', 'C': 'o', 'T': '^'}                     # hexagon / circle / triangle (richer shapes)
CMAP  = {'H': plt.cm.Oranges, 'C': plt.cm.Greens, 'T': plt.cm.Blues}   # within-class degree gradient
COL   = {'H': '#E69F00', 'C': '#2CA25F', 'T': '#3182BD'}   # legend representative mid-tone
SZ    = {'H': 1500, 'C': 150, 'T': 320}                    # base size (scaled by degree)
R     = {'H': 0.70, 'C': 1.70, 'T': 2.80}                  # concentric radii: herb(inner)/compound(mid)/target(outer)
groups = {k: [n for n in G if typ(n) == k] for k in 'HCT'}

def circmean(a):
    if not a: return 0.0
    return math.atan2(sum(math.sin(x) for x in a), sum(math.cos(x) for x in a))

# concentric angular placement, ordered so neighbours sit close (minimise crossings)
ang = {}
horder = [h for h in ['H:sanqi','H:chuanxiong','H:yanhusuo','H:qianghuo','H:baishao','H:weilingxian','H:gegen'] if h in groups['H']]
for i, h in enumerate(horder): ang[h] = 2*math.pi*i/len(horder)
comp = sorted(groups['C'], key=lambda c: circmean([ang[nb] for nb in G.neighbors(c) if nb in ang]) % (2*math.pi))
for i, c in enumerate(comp): ang[c] = 2*math.pi*i/len(comp)
tgt = sorted(groups['T'], key=lambda t: circmean([ang[nb] for nb in G.neighbors(t) if typ(nb)=='C']) % (2*math.pi))
for i, t in enumerate(tgt): ang[t] = 2*math.pi*i/len(tgt)
pos = {n: (R[typ(n)]*math.cos(ang[n]), R[typ(n)]*math.sin(ang[n])) for n in G}

fig, ax = plt.subplots(figsize=(8.4, 8.4)); ax.set_aspect('equal')
# faint concentric guide rings
for r in (R['C'], R['T']):
    ax.add_patch(plt.Circle((0, 0), r, fill=False, color='#E4E4E4', lw=0.8, zorder=0))
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.10, width=0.3, edge_color='#9AA0A6')
for k in 'HCT':                                            # shape per class + within-class degree gradient
    nodes = groups[k]
    dv = np.array([deg[n] for n in nodes], float)
    t = (dv - dv.min())/(dv.max()-dv.min()+1e-9) if len(dv) and dv.max() > dv.min() else np.full(len(dv), 0.6)
    colors = [CMAP[k](0.38 + 0.57*ti) for ti in t]
    sizes  = [SZ[k]*(0.5 + 1.2*ti) for ti in t]
    nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=colors, node_shape=SHAPE[k],
                           node_size=sizes, edgecolors='black', linewidths=0.4, ax=ax)
# herb labels (centre, white text with dark stroke -> readable on any orange shade)
for n in groups['H']:
    tx = ax.text(*pos[n], HERB_CN.get(n[2:], n[2:]), fontsize=7, ha='center', va='center',
                 fontweight='bold', color='white', zorder=5); F.stroke(tx)
# target labels (outer, radial; hubs bold)
for n in groups['T']:
    a = ang[n]; r = R['T'] + 0.13; da = math.degrees(a)
    rot = da if -90 < da < 90 else da + 180; ha = 'left' if -90 < da < 90 else 'right'
    g = n[2:]
    ax.text(r*math.cos(a), r*math.sin(a), g, fontsize=5.2, rotation=rot, rotation_mode='anchor',
            ha=ha, va='center', fontweight='bold' if g in hubs else 'normal',
            color='#000000' if g in hubs else '#666666')
# key compound labels (top degree, radial, inside the middle ring)
for n in sorted(groups['C'], key=lambda x: -deg[x])[:10]:
    a = ang[n]; r = R['C'] - 0.16; da = math.degrees(a)
    rot = da if -90 < da < 90 else da + 180; ha = 'right' if -90 < da < 90 else 'left'
    ax.text(r*math.cos(a), r*math.sin(a), n[2:][:18], fontsize=4.6, rotation=rot,
            rotation_mode='anchor', ha=ha, va='center', color='#0B3D2E')
ax.legend(handles=[Line2D([0],[0],marker=SHAPE[k],color='none',markerfacecolor=COL[k],markeredgecolor='black',
                          markersize=ms, label=lab)
                   for k, ms, lab in (('H',12,'Herb (inner, hexagon)'),('C',7,'Active compound (middle, circle)'),
                                      ('T',8,'Drug–disease target (outer, triangle)'))],
          loc='upper left', bbox_to_anchor=(-0.02, 1.0), fontsize=7.5, frameon=False)
ax.set_title('Herb–compound–target network of Jingtong Granules (7 herbs)', fontsize=10)
lim = R['T'] + 0.6; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.axis('off')
ax.text(0, -lim+0.12, 'Shape = node class; colour depth & size scale with degree within each class',
        ha='center', fontsize=6.6, color='#555555')
fig.tight_layout()
F.save(fig, os.path.join(FIG, 'Fig_herb_compound_target'))
# export for Cytoscape
pd.DataFrame([(a.split(':',1)[1], b.split(':',1)[1], typ(a)+'-'+typ(b)) for a,b in G.edges()],
             columns=['source','target','edge_type']).to_csv(os.path.join(OUT,'herb_compound_target_edges.csv'),index=False,encoding='utf-8-sig')
print('saved Fig_herb_compound_target + edges.csv | nodes', G.number_of_nodes(), 'edges', G.number_of_edges())
