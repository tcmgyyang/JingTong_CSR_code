# -*- coding: utf-8 -*-
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

# ---------------- shared style ----------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

# Okabe-Ito palette
BLUE = "#0072B2"; ORANGE = "#E69F00"; GREEN = "#009E73"
VERM = "#D55E00"; PURPLE = "#CC79A7"; SKY = "#56B4E9"
YELLOW = "#F0E442"; GREY = "#999999"; BLACK = "#000000"

BASE = r"H:/毕业设计/网药部分/JingTong_CSR_paper"
NODES_CSV = BASE + "/tables/cerna/ceRNA_nodes.csv"
EDGES_CSV = BASE + "/tables/cerna/ceRNA_edges.csv"
OUT = BASE + "/figures/composite/FigS3_ceRNA"

# ---------------- load ----------------
nodes = pd.read_csv(NODES_CSV)
edges = pd.read_csv(EDGES_CSV)

ntype = dict(zip(nodes["node"], nodes["type"]))

mrna  = [n for n in nodes["node"] if ntype[n] == "mRNA"]
mirna = [n for n in nodes["node"] if ntype[n] == "miRNA"]
lncrna = [n for n in nodes["node"] if ntype[n] == "lncRNA"]

# Fixed hub order requested
hub_order = ["IL1B", "IL6", "JUN", "BCL2L1", "CCL2", "CXCL8"]
mrna = [m for m in hub_order if m in mrna] + [m for m in mrna if m not in hub_order]

n_m, n_mi, n_lnc = len(mrna), len(mirna), len(lncrna)
print("counts mRNA/miRNA/lncRNA:", n_m, n_mi, n_lnc)

# ---------------- graph ----------------
G = nx.Graph()
for n in nodes["node"]:
    G.add_node(n, type=ntype[n])
for _, r in edges.iterrows():
    G.add_edge(r["source"], r["target"])
print("nodes/edges:", G.number_of_nodes(), G.number_of_edges())

# ---------------- concentric layout ----------------
# Order each ring to reduce edge crossings: place miRNAs near their mRNA hub,
# place lncRNAs near the angular mean of the miRNAs they connect to.

def ring_positions(items, radius, start=np.pi/2):
    pos = {}
    nN = len(items)
    for i, it in enumerate(items):
        ang = start + 2*np.pi * i / nN
        pos[it] = (radius*np.cos(ang), radius*np.sin(ang))
    return pos, {it: (start + 2*np.pi*i/nN) for i, it in enumerate(items)}

R_M, R_MI, R_LNC = 0.0, 3.2, 6.2  # central hubs clustered tight, rings outward

# place hubs in a small inner cluster (hexagon)
pos = {}
hub_r = 1.0
for i, m in enumerate(mrna):
    ang = np.pi/2 + 2*np.pi * i / n_m
    pos[m] = (hub_r*np.cos(ang), hub_r*np.sin(ang))
hub_ang = {m: np.pi/2 + 2*np.pi*i/n_m for i, m in enumerate(mrna)}

# Order miRNAs by the angle of the hub(s) they connect to -> cluster around hub
mi_anchor = {}
for mi in mirna:
    hubs = [nb for nb in G.neighbors(mi) if ntype.get(nb) == "mRNA"]
    if hubs:
        # circular mean of hub angles
        angs = [hub_ang[h] for h in hubs]
        s = np.mean([np.sin(a) for a in angs]); c = np.mean([np.cos(a) for a in angs])
        mi_anchor[mi] = np.arctan2(s, c)
    else:
        mi_anchor[mi] = 0.0
mirna_sorted = sorted(mirna, key=lambda x: mi_anchor[x])

mi_pos, mi_ang = ring_positions(mirna_sorted, R_MI)
pos.update(mi_pos)

# Order lncRNAs by angular mean of connected miRNAs
lnc_anchor = {}
for ln in lncrna:
    mis = [nb for nb in G.neighbors(ln) if ntype.get(nb) == "miRNA"]
    if mis:
        angs = [mi_ang[m] for m in mis if m in mi_ang]
        if angs:
            s = np.mean([np.sin(a) for a in angs]); c = np.mean([np.cos(a) for a in angs])
            lnc_anchor[ln] = np.arctan2(s, c)
        else:
            lnc_anchor[ln] = 0.0
    else:
        lnc_anchor[ln] = 0.0
lncrna_sorted = sorted(lncrna, key=lambda x: lnc_anchor[x])

lnc_pos, lnc_ang = ring_positions(lncrna_sorted, R_LNC)
pos.update(lnc_pos)

# ---------------- draw ----------------
fig, ax = plt.subplots(figsize=(9, 9))

# edges: thin grey
nx.draw_networkx_edges(
    G, pos, ax=ax,
    edge_color=GREY, width=0.45, alpha=0.5,
)

# node draws by class
nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=lncrna,
                       node_color=SKY, node_size=70,
                       edgecolors="white", linewidths=0.4)
nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=mirna,
                       node_color=GREEN, node_size=130,
                       edgecolors="white", linewidths=0.5)
nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=mrna,
                       node_color=VERM, node_size=620,
                       edgecolors="white", linewidths=1.2)

# ----- labels -----
# mRNA hubs: always labelled, bold
for m in mrna:
    x, y = pos[m]
    ax.text(x, y, m, fontsize=8.5, fontweight="bold", ha="center", va="center",
            color="white", zorder=6)

# miRNA labels: small, rotated radially, stripped 'hsa-' prefix
for mi in mirna_sorted:
    x, y = pos[mi]
    ang = np.degrees(mi_ang[mi])
    rot = ang
    ha = "left"
    if 90 < ang % 360 < 270:
        rot = ang + 180
        ha = "right"
    lbl = mi.replace("hsa-miR-", "miR-").replace("hsa-let-", "let-")
    # push label slightly outside the miRNA ring
    rx, ry = (R_MI + 0.28) * np.cos(mi_ang[mi]), (R_MI + 0.28) * np.sin(mi_ang[mi])
    ax.text(rx, ry, lbl, fontsize=4.6, ha=ha, va="center",
            rotation=rot, rotation_mode="anchor", color="#0b3d2e", zorder=5)

# lncRNA labels: outer ring, small grey-blue text radially
for ln in lncrna_sorted:
    ang = np.degrees(lnc_ang[ln])
    rot = ang
    ha = "left"
    if 90 < ang % 360 < 270:
        rot = ang + 180
        ha = "right"
    rx, ry = (R_LNC + 0.22) * np.cos(lnc_ang[ln]), (R_LNC + 0.22) * np.sin(lnc_ang[ln])
    ax.text(rx, ry, ln, fontsize=4.2, ha=ha, va="center",
            rotation=rot, rotation_mode="anchor", color="#15506e", zorder=5)

# ---------------- legend ----------------
legend_elems = [
    Line2D([0], [0], marker="o", linestyle="none", markersize=13,
           markerfacecolor=VERM, markeredgecolor="white",
           label=f"mRNA hub (n={n_m})"),
    Line2D([0], [0], marker="o", linestyle="none", markersize=8,
           markerfacecolor=GREEN, markeredgecolor="white",
           label=f"miRNA (n={n_mi})"),
    Line2D([0], [0], marker="o", linestyle="none", markersize=6,
           markerfacecolor=SKY, markeredgecolor="white",
           label=f"lncRNA (n={n_lnc})"),
    Line2D([0], [0], color=GREY, lw=0.8, alpha=0.7,
           label=f"ceRNA interaction (n={G.number_of_edges()})"),
]
ax.legend(handles=legend_elems, loc="upper left", bbox_to_anchor=(-0.02, 1.0),
          fontsize=7.5, handletextpad=0.6, labelspacing=0.9, borderpad=0.4)

ax.set_title("ceRNA network (lncRNA-miRNA-mRNA)", fontsize=13, fontweight="bold", pad=12)

ax.set_aspect("equal")
ax.axis("off")
lim = R_LNC + 1.7
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

fig.tight_layout()

for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")

# report pixel size of PNG
try:
    from PIL import Image
    w, h = Image.open(f"{OUT}.png").size
    print("PNG pixel size:", w, "x", h)
except Exception as e:
    print("size check failed:", e)

plt.close(fig)
print("DONE")
