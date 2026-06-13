# Supplementary figures — native: S1 immune, S2 context, S3 ceRNA network, S4 ADMET.
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"; os.environ["OMP_NUM_THREADS"] = "1"
import csv, traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import Patch, Ellipse
from matplotlib.lines import Line2D
try:
    from adjustText import adjust_text
except Exception:
    adjust_text = None

ROOT = r"H:\毕业设计\网药部分\JingTong_CSR_paper"
ML = os.path.join(ROOT, "tables", "ML"); CER = os.path.join(ROOT, "tables", "cerna")
ADM = os.path.join(ROOT, "tables", "admet"); CMP = os.path.join(ROOT, "figures", "composite")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.8, "svg.fonttype": "none", "pdf.fonttype": 42})
OK = {"CON": "#56B4E9", "CSR": "#E69F00", "DCM": "#D55E00"}
BLUE, VERM, GREEN, SKY, PURPLE, GREY = "#0072B2", "#D55E00", "#009E73", "#56B4E9", "#CC79A7", "#cccccc"
PAL12 = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442",
         "#999999", "#882255", "#44AA99", "#332288", "#AA4499"]
lab = pd.read_csv(os.path.join(ML, "labels.csv"), index_col=0).iloc[:, 0]
ORD = ["CON", "CSR", "DCM"]


def L(ax, t, x=-0.13, y=1.05):
    ax.text(x, y, t, transform=ax.transAxes, fontsize=15, fontweight="bold", va="bottom")


def save(fig, name):
    fig.savefig(os.path.join(CMP, f"{name}.png"), dpi=200, bbox_inches="tight", facecolor="white")
    for ext in ["pdf", "svg"]:
        try:
            fig.savefig(os.path.join(CMP, f"{name}.{ext}"), bbox_inches="tight", facecolor="white")
        except Exception:
            pass
    plt.close(fig); print("wrote", name)


def cat_grad(nodes, deg, cmap_name, lo=0.3):
    cmap = matplotlib.colormaps[cmap_name]
    vals = [deg[n] for n in nodes]
    lo_v, hi_v = (min(vals), max(vals)) if vals else (0, 1)
    span = (hi_v - lo_v) or 1
    return [cmap(lo + (1 - lo) * ((deg[n] - lo_v) / span)) for n in nodes]


def rlabel(ax, name, x, y, off, fs, color):
    ang = np.degrees(np.arctan2(y, x)); r = np.hypot(x, y) + off
    rot = ang if -90 <= ang <= 90 else ang + 180
    ha = "left" if -90 <= ang <= 90 else "right"
    ax.text(r * np.cos(np.radians(ang)), r * np.sin(np.radians(ang)), name, fontsize=fs,
            rotation=rot, rotation_mode="anchor", ha=ha, va="center", color=color)


# ---------------- S1 immune ----------------
def s1():
    cib = pd.read_csv(os.path.join(ML, "CIBERSORT_fractions.csv"), index_col=0)
    cib = cib.loc[[s for s in lab.index if s in cib.index]]
    grp = lab.loc[cib.index]
    LIN = {"B cells naive": "B", "B cells memory": "B", "Plasma cells": "Plasma",
           "T cells CD8": "CD8 T", "T cells CD4 naive": "CD4 T", "T cells CD4 memory resting": "CD4 T",
           "T cells CD4 memory activated": "CD4 T", "T cells follicular helper": "CD4 T",
           "T cells regulatory (Tregs)": "Treg", "T cells gamma delta": "γδ T",
           "NK cells resting": "NK", "NK cells activated": "NK", "Monocytes": "Monocyte",
           "Macrophages M0": "Macrophage", "Macrophages M1": "Macrophage", "Macrophages M2": "Macrophage",
           "Dendritic cells resting": "DC", "Dendritic cells activated": "DC",
           "Mast cells resting": "Mast", "Mast cells activated": "Mast",
           "Eosinophils": "Eosinophil", "Neutrophils": "Neutrophil"}
    coll = cib.rename(columns=LIN).groupby(level=0, axis=1).sum()
    order = grp.sort_values().index
    lineages = ["Neutrophil", "CD8 T", "CD4 T", "Monocyte", "Macrophage", "NK", "γδ T",
                "B", "Treg", "Mast", "DC", "Plasma", "Eosinophil"]
    lineages = [x for x in lineages if x in coll.columns]

    ss = pd.read_csv(os.path.join(ML, "ssGSEA_immune.csv"), index_col=0)
    ss = ss.loc[[s for s in lab.index if s in ss.index]]
    sst = {r["cell"]: float(r["p_CERVvsCON"]) for r in csv.DictReader(open(os.path.join(ML, "ssGSEA_group_stats.csv"), encoding="utf-8-sig"))}

    fig = plt.figure(figsize=(12.5, 9.2)); gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.26)
    # A stacked bar
    axA = fig.add_subplot(gs[0, 0]); bottom = np.zeros(len(order))
    for i, ln in enumerate(lineages):
        axA.bar(range(len(order)), coll.loc[order, ln].values, bottom=bottom,
                color=PAL12[i % len(PAL12)], width=1.0, label=ln)
        bottom += coll.loc[order, ln].values
    axA.set_xlim(-0.5, len(order) - 0.5); axA.set_ylabel("Fraction", fontsize=8.5)
    axA.set_title("CIBERSORT immune composition", fontsize=9.5)
    axA.set_xticks([]); axA.legend(fontsize=6, ncol=2, loc="upper right", bbox_to_anchor=(1.0, 1.0))
    # group boundaries
    gc = grp.loc[order].values
    for k in ORD:
        idx = np.where(gc == k)[0]
        if len(idx): axA.text(idx.mean(), -0.05, k, ha="center", va="top", fontsize=7, transform=axA.get_xaxis_transform())
    L(axA, "A")
    # B ssGSEA significant boxplots
    axB = fig.add_subplot(gs[0, 1]); sig = ["CD8_T_cell", "NK_cell", "Cytotoxic", "Mast_cell"]
    sig = [c for c in sig if c in ss.columns]; w = 0.24
    for gi, c in enumerate(sig):
        for ki, k in enumerate(ORD):
            d = ss.loc[grp[grp == k].index, c].values
            bp = axB.boxplot([d], positions=[gi + (ki - 1) * w], widths=w * 0.85, patch_artist=True,
                             showfliers=False, manage_ticks=False)
            bp["boxes"][0].set(facecolor=OK[k], alpha=0.85, edgecolor="black", linewidth=0.5)
            bp["medians"][0].set(color="black", linewidth=0.8)
        axB.text(gi, axB.get_ylim()[1], f"p={sst.get(c,float('nan')):.3f}", ha="center", va="bottom", fontsize=6.5)
    axB.set_xticks(range(len(sig))); axB.set_xticklabels([c.replace("_", " ") for c in sig], fontsize=7.5)
    axB.set_ylabel("ssGSEA score", fontsize=8.5); axB.set_title("ssGSEA: altered immune programs", fontsize=9.5)
    axB.legend(handles=[Patch(facecolor=OK[k], label=k) for k in ORD], fontsize=7, loc="lower right")
    axB.spines[["top", "right"]].set_visible(False); L(axB, "B")
    # C -log10 p
    axC = fig.add_subplot(gs[1, 0]); cells = sorted(sst, key=lambda c: sst[c])
    vals = [-np.log10(sst[c]) for c in cells]
    axC.barh(range(len(cells)), vals, color=[VERM if sst[c] < 0.05 else GREY for c in cells])
    axC.axvline(-np.log10(0.05), color="0.5", ls="--", lw=0.8)
    axC.set_yticks(range(len(cells))); axC.set_yticklabels([c.replace("_", " ") for c in cells], fontsize=7.5)
    axC.invert_yaxis(); axC.set_xlabel("−log10 p (degeneration vs control)", fontsize=8.5)
    axC.set_title("ssGSEA group differences", fontsize=9.5); axC.spines[["top", "right"]].set_visible(False); L(axC, "C")
    # D CIBERSORT top cells boxplots
    axD = fig.add_subplot(gs[1, 1]); top = ["Neutrophil", "CD8 T", "NK", "Monocyte"]
    top = [t for t in top if t in coll.columns]
    for gi, c in enumerate(top):
        for ki, k in enumerate(ORD):
            d = coll.loc[grp[grp == k].index, c].values
            bp = axD.boxplot([d], positions=[gi + (ki - 1) * w], widths=w * 0.85, patch_artist=True,
                             showfliers=False, manage_ticks=False)
            bp["boxes"][0].set(facecolor=OK[k], alpha=0.85, edgecolor="black", linewidth=0.5)
            bp["medians"][0].set(color="black", linewidth=0.8)
    axD.set_xticks(range(len(top))); axD.set_xticklabels(top, fontsize=7.5)
    axD.set_ylabel("CIBERSORT fraction", fontsize=8.5); axD.set_title("CIBERSORT cell fractions", fontsize=9.5)
    axD.spines[["top", "right"]].set_visible(False); L(axD, "D")
    save(fig, "FigS1_immune")


# ---------------- S2 context ----------------
def s2():
    tis = pd.read_csv(os.path.join(ML, "tissue_hubgene_GSE153761.csv"))
    gcol = tis.columns[0]
    fccol = [c for c in tis.columns if "log2" in c.lower() or "fc" in c.lower()][0]
    pcol = [c for c in tis.columns if c.lower() in ("p", "pval", "pvalue") or c.lower().startswith("p_")][0]
    tis = tis.sort_values(fccol)
    meta = pd.read_csv(os.path.join(ML, "sample_meta.csv"))
    age = pd.read_csv(os.path.join(ML, "hubgene_age_correlation.csv"))
    himm = pd.read_csv(os.path.join(ML, "hub_immune_spearman_r.csv"), index_col=0)

    fig = plt.figure(figsize=(12.5, 9.0)); gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.30)
    # A tissue log2FC
    axA = fig.add_subplot(gs[0, 0]); fc = tis[fccol].values
    axA.barh(range(len(tis)), fc, color=[VERM if v >= 0 else BLUE for v in fc], edgecolor="black", linewidth=0.4)
    for i, (v, p) in enumerate(zip(fc, tis[pcol].values)):
        if p < 0.05: axA.text(v + (0.03 if v >= 0 else -0.03), i, "*", ha="center", va="center", fontsize=12, color="black")
    axA.set_yticks(range(len(tis))); axA.set_yticklabels(tis[gcol].values, fontsize=7.5)
    axA.axvline(0, color="0.6", lw=0.8); axA.set_xlabel("log2 FC (degenerate vs normal)", fontsize=8.5)
    axA.set_title("Tissue (GSE153761): * p<0.05", fontsize=9.5); axA.spines[["top", "right"]].set_visible(False); L(axA, "A")
    # B age confound: left = age by group (the confound); right = hub~age correlation (not the driver)
    subB = gs[0, 1].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.62)
    axB = fig.add_subplot(subB[0]); axB2 = fig.add_subplot(subB[1])
    gcol2 = "group" if "group" in meta.columns else meta.columns[1]
    acol = "age" if "age" in meta.columns else [c for c in meta.columns if "age" in c.lower()][0]
    for ki, k in enumerate(ORD):
        d = meta[meta[gcol2] == k][acol].values
        axB.scatter(np.random.RandomState(ki).normal(ki, 0.06, len(d)), d, color=OK[k], s=18, alpha=0.8, zorder=3)
        bp = axB.boxplot([d], positions=[ki], widths=0.5, patch_artist=True, showfliers=False, manage_ticks=False)
        bp["boxes"][0].set(facecolor=OK[k], alpha=0.3, edgecolor=OK[k]); bp["medians"][0].set(color="black")
    axB.set_xticks(range(3)); axB.set_xticklabels(ORD, fontsize=8); axB.set_ylabel("Age (years)", fontsize=8.5)
    axB.set_title("Controls younger", fontsize=9); axB.spines[["top", "right"]].set_visible(False); L(axB, "B", x=-0.32)
    # right: hub-gene ~ age correlation (all weak -> age is not driving the signal)
    rcol = [c for c in age.columns if "spearman" in c.lower() or c.lower() == "r"][0]
    ag = age.sort_values(rcol); yy = np.arange(len(ag))
    axB2.axvspan(-0.2, 0.2, color="#ededed", zorder=0)
    axB2.hlines(yy, 0, ag[rcol].values, color="#9aa0a6", lw=2, zorder=2)
    axB2.scatter(ag[rcol].values, yy, color=PURPLE, s=40, zorder=3)
    axB2.axvline(0, color="0.55", lw=0.8)
    axB2.set_yticks(yy); axB2.set_yticklabels(ag[age.columns[0]].values, fontsize=7.5)
    axB2.set_xlim(-0.55, 0.55); axB2.set_xlabel("Spearman r vs age", fontsize=8)
    axB2.set_title("Weak (|r|≤0.2):\nage not the driver", fontsize=8); axB2.spines[["top", "right"]].set_visible(False)
    # C hub-immune correlation heatmap
    axC = fig.add_subplot(gs[1, :])
    M = himm.values.astype(float)
    im = axC.imshow(M, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    axC.set_yticks(range(himm.shape[0])); axC.set_yticklabels(himm.index, fontsize=8)
    axC.set_xticks(range(himm.shape[1])); axC.set_xticklabels(himm.columns, rotation=55, ha="right", fontsize=6.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isnan(v):
                axC.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.0,
                         color="white" if abs(v) > 0.35 else "black")
    cb = fig.colorbar(im, ax=axC, fraction=0.02, pad=0.01); cb.set_label("Spearman r", fontsize=8)
    axC.set_title("Hub-gene – immune-cell correlation", fontsize=9.5); L(axC, "C", x=-0.06)
    save(fig, "FigS2_context")


# ---------------- S3 ceRNA concentric ----------------
def s3():
    ntype = {r["node"]: r["type"] for r in csv.DictReader(open(os.path.join(CER, "ceRNA_nodes.csv"), encoding="utf-8-sig"))}
    G = nx.Graph()
    for r in csv.DictReader(open(os.path.join(CER, "ceRNA_edges.csv"), encoding="utf-8-sig")):
        G.add_edge(r["source"], r["target"])
    mR = [n for n in ntype if ntype[n] == "mRNA"]; miR = [n for n in ntype if ntype[n] == "miRNA"]
    lnc = [n for n in ntype if ntype[n] == "lncRNA"]

    def ring(nodes, rad, start=90.0):
        return {nd: (rad * np.cos(np.deg2rad(start - 360 * i / max(len(nodes), 1))),
                     rad * np.sin(np.deg2rad(start - 360 * i / max(len(nodes), 1)))) for i, nd in enumerate(nodes)}
    pos = {}; pos.update(ring(mR, 0.28)); pos.update(ring(miR, 0.62)); pos.update(ring(lnc, 1.0))
    fig, ax = plt.subplots(figsize=(11.5, 11.5))
    deg = dict(G.degree())
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#ececec", width=0.32)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=lnc, node_color=cat_grad(lnc, deg, "Blues", 0.30),
                           node_size=[22 + deg[n] * 9 for n in lnc], edgecolors="white", linewidths=0.3)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=miR, node_color=cat_grad(miR, deg, "Greens", 0.35),
                           node_size=[30 + deg[n] * 6 for n in miR], edgecolors="white", linewidths=0.3)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=mR, node_color=cat_grad(mR, deg, "OrRd", 0.5),
                           node_size=560, edgecolors="black", linewidths=0.7)

    def shorten(n, k=13):
        n = n.replace("hsa-", "")
        return n if len(n) <= k else n[:k - 1] + "…"
    for n in lnc:
        rlabel(ax, shorten(n), pos[n][0], pos[n][1], 0.02, 3.8, "#1f5fa6")
    for n in miR:
        rlabel(ax, shorten(n), pos[n][0], pos[n][1], 0.015, 3.4, "#176b4a")
    for n in mR:
        ax.text(*pos[n], n, fontsize=7.5, ha="center", va="center", fontweight="bold", color="white")
    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=VERM, markersize=11, label=f"mRNA hub ({len(mR)}, inner)"),
                       Line2D([0], [0], marker="o", color="w", markerfacecolor=GREEN, markersize=8, label=f"miRNA ({len(miR)}, middle)"),
                       Line2D([0], [0], marker="o", color="w", markerfacecolor=SKY, markersize=7, label=f"lncRNA ({len(lnc)}, outer)")],
              fontsize=8, loc="upper right", frameon=False)
    ax.set_xlim(-1.42, 1.42); ax.set_ylim(-1.42, 1.42); ax.set_aspect("equal"); ax.axis("off")
    ax.margins(0)
    ax.set_title("ceRNA network (lncRNA–miRNA–mRNA)", fontsize=11)
    save(fig, "FigS3_ceRNA")


# ---------------- S4 ADMET ----------------
def s4():
    df = pd.read_csv(os.path.join(ADM, "swissadme.csv"))
    name = df.columns[0]
    props = {"MW": "MW", "Consensus Log P": "cLogP", "TPSA": "TPSA", "#Rotatable bonds": "RotB",
             "#H-bond donors": "HBD", "#H-bond acceptors": "HBA", "Synthetic Accessibility": "SA"}
    props = {k: v for k, v in props.items() if k in df.columns}
    M = df[list(props)].astype(float).values
    Mn = (M - M.min(0)) / (M.max(0) - M.min(0) + 1e-9)
    fig = plt.figure(figsize=(13.5, 5.8)); gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.28)
    axA = fig.add_subplot(gs[0]); im = axA.imshow(Mn, cmap="RdBu_r", vmin=0, vmax=1, aspect="auto")
    axA.set_xticks(range(len(props))); axA.set_xticklabels(list(props.values()), fontsize=8)
    axA.set_yticks(range(len(df))); axA.set_yticklabels(df[name].values, fontsize=8)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            axA.text(j, i, f"{M[i,j]:g}", ha="center", va="center", fontsize=6,
                     color="white" if abs(Mn[i, j] - 0.5) > 0.32 else "black")
    cb = fig.colorbar(im, ax=axA, fraction=0.045, pad=0.02); cb.set_label("min–max normalised", fontsize=8)
    axA.set_title("Physicochemical properties (SwissADME)", fontsize=9.5); L(axA, "A", x=-0.18)
    # B BOILED-Egg
    axB = fig.add_subplot(gs[1])
    tpsa = df["TPSA"].astype(float).values; wlogp = df["WLOGP"].astype(float).values
    axB.add_patch(Ellipse((71.05, 2.292), 2 * 142.08, 2 * 8.74, color="#f3f3f3", zorder=0))
    axB.add_patch(Ellipse((38.117, 3.177), 2 * 82.06, 2 * 5.557, color="#fde9a9", zorder=1))
    gi = df["GI absorption"].astype(str).values if "GI absorption" in df.columns else np.array(["?"] * len(df))
    for i in range(len(df)):
        c = VERM if gi[i] == "High" else BLUE
        axB.scatter(tpsa[i], wlogp[i], color=c, s=44, edgecolors="black", linewidths=0.4, zorder=3)
    # widen the view so labels have room, then repel them off each other and the points (leader lines)
    axB.set_xlim(min(tpsa) - 60, max(tpsa) + 60); axB.set_ylim(min(wlogp) - 2.4, max(wlogp) + 2.6)
    txts = [axB.text(tpsa[i], wlogp[i], df[name].values[i], fontsize=6, zorder=5) for i in range(len(df))]
    if adjust_text is not None:
        adjust_text(txts, x=list(tpsa), y=list(wlogp), ax=axB,
                    expand=(1.25, 1.6), force_text=(0.5, 0.9), force_static=(0.3, 0.5),
                    arrowprops=dict(arrowstyle="-", color="0.45", lw=0.5, shrinkA=3, shrinkB=4))
    else:  # fallback: fixed offset (overlap-prone)
        for t, (tx, ty) in zip(txts, zip(tpsa, wlogp)):
            t.set_position((tx + 2, ty + 2)); t.set_fontsize(5.5)
    axB.set_xlabel("TPSA (Å²)", fontsize=8.5); axB.set_ylabel("WLOGP", fontsize=8.5)
    axB.set_title("BOILED-Egg (yolk = BBB, white = HIA)", fontsize=9.5)
    axB.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=VERM, markersize=8, label="GI High"),
                        Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE, markersize=8, label="GI Low")],
               fontsize=7.5, loc="upper right"); axB.spines[["top", "right"]].set_visible(False); L(axB, "B", x=-0.14)
    save(fig, "FigS4_admet")


for fn in [s1, s2, s3, s4]:
    try:
        fn()
    except Exception:
        print("FAILED", fn.__name__); traceback.print_exc()
    finally:
        plt.close("all")
