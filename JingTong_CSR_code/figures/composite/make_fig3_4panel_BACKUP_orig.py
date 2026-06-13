# -*- coding: utf-8 -*-
"""Fig 3 (4-panel): a = GO-BP ORA, b = KEGG ORA, c = GSVA, d = single-gene GSEA.
GO/KEGG come from Enrichr ORA on the 43 intersection targets (tables/enrichment/*_ORA.csv);
GSVA/GSEA reuse the existing key-pathway summaries. Reproduces Results 3.3."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
    "font.size":8,"axes.linewidth":0.8,"axes.spines.top":False,"axes.spines.right":False,
    "legend.frameon":False,"svg.fonttype":"none","pdf.fonttype":42})
BLUE="#0072B2"; ORANGE="#E69F00"; GREEN="#009E73"; VERM="#D55E00"; GREY="#999999"; BLACK="#000000"

BASE = r"O:\JingTong_CSR_paper"
ENR = os.path.join(BASE,"tables","enrichment")
OUT = os.path.join(BASE,"figures","composite","Fig3_pathways_4panel")

keg = pd.read_csv(os.path.join(ENR,"KEGG_ORA.csv"))
gsva= pd.read_csv(os.path.join(ENR,"GSVA_key_pathways.csv"))
gsea= pd.read_csv(os.path.join(ENR,"GSEA_key_pathways_summary.csv"))
from matplotlib.patches import Patch
SKY="#56B4E9"
GO_FILES={"BP":"GO_BP_ORA.csv","CC":"GO_CC_ORA.csv","MF":"GO_MF_ORA.csv"}
GO_COLORS={"BP":GREEN,"CC":SKY,"MF":ORANGE}
GO_NAMES={"BP":"Biological process","CC":"Cellular component","MF":"Molecular function"}

# KEGG: drop promiscuous disease/infection/cancer gene-sets -> keep signaling/mechanism pathways
BAD = ("infection","disease","cancer","sarcoma","leukemia","carcinoma","hepatitis",
       "influenza","tuberculosis","measles","Epstein","virus","Yersinia","Chagas",
       "Salmonella","Shigellosis","Amoebiasis","Toxoplasmosis","Malaria","Pertussis",
       "Legionellosis","Leishmaniasis","platelet")
keg_sig = keg[~keg["term"].str.lower().apply(lambda t: any(b.lower() in t for b in BAD))].copy()

def clean_go(t): return t.rsplit(" (GO:",1)[0]

fig = plt.figure(figsize=(12.6, 10.4))
gs = GridSpec(2, 2, figure=fig, wspace=0.62, hspace=0.42,
              left=0.085, right=0.93, top=0.93, bottom=0.08,
              height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.05])
axA = fig.add_subplot(gs[0,0]); axB = fig.add_subplot(gs[0,1])
axC = fig.add_subplot(gs[1,0]); axD = fig.add_subplot(gs[1,1])

def hbar(ax, df, termcol, color, title, panel, n=10, label_clean=None):
    d = df.sort_values("adj_p").head(n).copy()
    d["nlp"] = -np.log10(d["adj_p"].clip(lower=1e-300))
    terms = d[termcol].tolist()
    if label_clean: terms = [label_clean(t) for t in terms]
    y = np.arange(len(d))[::-1]  # most-significant at top
    ax.barh(y, d["nlp"].values, color=color, edgecolor="white", linewidth=0.6, height=0.72, zorder=3)
    # gene-count annotation at bar end
    for yi, nlp, ncol in zip(y, d["nlp"].values, d["n_overlap"].values):
        ax.text(nlp + d["nlp"].max()*0.015, yi, str(int(ncol)), va="center", ha="left",
                fontsize=6.5, color=GREY)
    ax.set_yticks(y); ax.set_yticklabels(terms, fontsize=7.3)
    ax.set_xlabel("$-\\log_{10}$ adjusted $P$", fontsize=8)
    ax.set_xlim(0, d["nlp"].max()*1.16)
    ax.tick_params(axis="both", length=3, width=0.8)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=8)
    ax.text(-0.02, 1.05, panel, transform=ax.transAxes, fontsize=15, fontweight="bold", ha="right", va="bottom")

# Panel A: combined GO (BP / CC / MF), top 4 each, color-coded by ontology
rows = []  # (ontology, label, nlp, n_overlap)
for onto in ("BP", "CC", "MF"):
    d = pd.read_csv(os.path.join(ENR, GO_FILES[onto])).sort_values("adj_p").head(4)
    for _, r in d.iterrows():
        rows.append((onto, clean_go(r["term"]),
                     -np.log10(max(r["adj_p"], 1e-300)), int(r["n_overlap"])))
# assign y top->bottom (BP top), with a gap between ontology blocks
cur, yv, last = 0.0, [], None
for onto, *_ in rows:
    if last is not None and onto != last: cur += 0.9
    yv.append(cur); cur += 1.0; last = onto
ymax = max(yv); yv = [ymax - v for v in yv]
vals = [r[2] for r in rows]; cols = [GO_COLORS[r[0]] for r in rows]; ncs = [r[3] for r in rows]
axA.barh(yv, vals, color=cols, edgecolor="white", linewidth=0.6, height=0.74, zorder=3)
for yi, v, nc in zip(yv, vals, ncs):
    axA.text(v + max(vals)*0.015, yi, str(nc), va="center", ha="left", fontsize=6.3, color=GREY)
axA.set_yticks(yv); axA.set_yticklabels([r[1] for r in rows], fontsize=6.9)
axA.set_xlim(0, max(vals)*1.16); axA.set_ylim(min(yv)-0.7, max(yv)+0.7)
axA.set_xlabel("$-\\log_{10}$ adjusted $P$", fontsize=8)
axA.tick_params(axis="both", length=3, width=0.8)
axA.set_title("GO enrichment", fontsize=10, fontweight="bold", pad=8)
axA.text(-0.02, 1.05, "A", transform=axA.transAxes, fontsize=15, fontweight="bold", ha="right", va="bottom")
axA.legend(handles=[Patch(facecolor=GO_COLORS[o], edgecolor="none", label=GO_NAMES[o]) for o in ("BP","CC","MF")],
           loc="lower right", fontsize=6.6, frameon=False, handlelength=1.0, handleheight=1.0, labelspacing=0.35)
# Panel B: KEGG (signaling)
hbar(axB, keg_sig, "term", VERM, "KEGG pathway", "B", n=10)

# Panel C: GSVA lollipop
gC = gsva.sort_values("delta_CERVvsCON", ascending=True).reset_index(drop=True)
yc = np.arange(len(gC)); delta = gC["delta_CERVvsCON"].values; pv = gC["p"].values
cols = [VERM if v>=0 else BLUE for v in delta]
for yi,di,ci in zip(yc,delta,cols): axC.plot([0,di],[yi,yi],color=ci,lw=2.2,zorder=2,solid_capstyle="round")
axC.scatter(delta,yc,color=cols,s=52,zorder=3,edgecolor="white",linewidth=0.6)
axC.axvline(0,color=BLACK,lw=0.8)
axC.set_yticks(yc); axC.set_yticklabels(gC["pathway"].values,fontsize=7.3)
axC.set_xlabel("GSVA Δ (degeneration − control)",fontsize=8)
xpad=max((delta.max()-delta.min())*0.04,0.006)
for yi,di,p in zip(yc,delta,pv):
    t=f"p={p:.3g}"
    axC.text(di+(xpad if di>=0 else -xpad),yi,t,va="center",ha=("left" if di>=0 else "right"),fontsize=6.3,color=GREY)
sp=delta.max()-delta.min()
axC.set_xlim(min(delta.min(),0)-sp*0.42,max(delta.max(),0)+sp*0.42); axC.set_ylim(-0.7,len(gC)-0.3)
axC.tick_params(axis="both",length=3,width=0.8)
axC.set_title("GSVA pathway activity",fontsize=10,fontweight="bold",pad=8)
axC.text(-0.02,1.05,"C",transform=axC.transAxes,fontsize=15,fontweight="bold",ha="right",va="bottom")

# Panel D: single-gene GSEA NES heatmap
nes=gsea.pivot(index="gene",columns="Term",values="NES"); fdr=gsea.pivot(index="gene",columns="Term",values="FDR")
go_order=nes.mean(1).sort_values(ascending=False).index; nes=nes.loc[go_order]; fdr=fdr.loc[go_order]
to=nes.mean(0).sort_values(ascending=False).index; nes=nes[to]; fdr=fdr[to]
data=nes.values.astype(float); fdrv=fdr.values.astype(float)
half=max(1.0-np.nanmin(data),np.nanmax(data)-1.0); norm=TwoSlopeNorm(vmin=1.0-half,vcenter=1.0,vmax=1.0+half)
nr,nc=data.shape
im=axD.pcolormesh(np.arange(nc+1),np.arange(nr+1),data,cmap="RdBu_r",norm=norm,edgecolors="white",linewidth=1.0,shading="flat")
axD.set_aspect("auto"); axD.set_xlim(0,nc); axD.set_ylim(0,nr); axD.invert_yaxis()
xc=np.arange(nc)+0.5; yc2=np.arange(nr)+0.5
axD.set_xticks(xc); axD.set_yticks(yc2)
axD.set_xticklabels([t.replace(" signaling pathway","").replace(" in diabetic complications","") for t in nes.columns],fontsize=6.8,rotation=40,ha="right",rotation_mode="anchor")
axD.set_yticklabels(nes.index,fontsize=7.3); axD.tick_params(length=0)
for s in axD.spines.values(): s.set_visible(False)
for i in range(nr):
    for j in range(nc):
        if not np.isnan(fdrv[i,j]) and fdrv[i,j]<0.05:
            axD.text(xc[j],yc2[i],"*",ha="center",va="center",fontsize=11,color=BLACK,fontweight="bold")
axD.set_title("Single-gene GSEA (NES)",fontsize=10,fontweight="bold",pad=8)
axD.text(-0.02,1.05,"D",transform=axD.transAxes,fontsize=15,fontweight="bold",ha="right",va="bottom")
cb=fig.colorbar(im,ax=axD,fraction=0.046,pad=0.02); cb.set_label("NES",fontsize=8)
cb.ax.tick_params(labelsize=7,length=2); cb.outline.set_linewidth(0.6); cb.ax.axhline(norm(1.0),color=BLACK,lw=0.8)

for ext in ("png","pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
from PIL import Image
with Image.open(f"{OUT}.png") as im_: print("PNG_SIZE",im_.size)
print("KEGG panel terms:", keg_sig.sort_values('adj_p').head(10)['term'].tolist())
print("DONE ->", OUT+".png")
