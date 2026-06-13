# -*- coding: utf-8 -*-
"""
Route B - Blood immune-infiltration analysis of GSE223227 for the Jingtong Granules / CSR paper.
- CIBERSORT (LM22) immune-cell deconvolution on peripheral-blood RNA-seq
- Group comparison (CSR / DCM / CON)
- ROC of blood-expressed hub genes (IL1B, MMP9, AKT1, TP53) for cervical degeneration vs control
- Spearman correlation of hub genes with immune-cell fractions
All figures saved at 300 dpi (PNG + PDF) in ./figures
"""
import os, sys, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.svm import NuSVR
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc
from scipy import stats

# shared publication style (figstyle.py at project root)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F

ROOT = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(ROOT, 'figures'); os.makedirs(FIG, exist_ok=True)
F.apply()                      # global rcParams: Arial, 8pt, despined, colourblind-safe
PAL = F.PAL                    # Okabe-Ito group triad {'CON','CSR','DCM'}

def save(fig, name):
    F.save(fig, os.path.join(FIG, name))   # PNG (300 dpi) + PDF (fonts embedded)

# ---------- load ----------
expr = pd.read_csv(os.path.join(ROOT, 'expr_symbol.csv.gz'), index_col=0)        # genes x samples, FPKM
labels = pd.read_csv(os.path.join(ROOT, 'labels.csv'), index_col=0)['group']
expr = expr[labels.index]
log2 = np.log2(expr + 1)
print('expr', expr.shape, '| groups', labels.value_counts().to_dict())

# ================= 1. CIBERSORT (LM22) =================
def cibersort(signature, mixture, nus=(0.25, 0.5, 0.75)):
    genes = signature.index.intersection(mixture.index)
    sig = signature.loc[genes].astype(float)
    mix = mixture.loc[genes].astype(float)
    sig_z = (sig - sig.values.mean()) / sig.values.std()
    out = {}
    for s in mix.columns:
        y = mix[s].values
        y_z = (y - y.mean()) / y.std()
        best = None
        for nu in nus:
            m = NuSVR(kernel='linear', nu=nu, C=1.0).fit(sig_z.values, y_z)
            w = m.coef_.flatten(); w[w < 0] = 0
            if w.sum() > 0: w = w / w.sum()
            rmse = np.sqrt(np.mean((sig_z.values.dot(w) - y_z) ** 2))
            if best is None or rmse < best[0]: best = (rmse, w)
        out[s] = best[1]
    res = pd.DataFrame(out, index=sig.columns).T
    return res, len(genes)

lm22 = pd.read_csv(os.path.join(ROOT, 'LM22.txt'), sep='\t', index_col=0)
frac, ngene = cibersort(lm22, expr)
frac.to_csv(os.path.join(ROOT, 'CIBERSORT_fractions.csv'))
print(f'CIBERSORT done: {ngene} shared LM22 genes; fractions {frac.shape}')

# ---- Fig 1: stacked bar of immune fractions ----
# 22 LM22 cells are non-separable by any qualitative colormap; collapse to ~10 lineages
# by SUMMING fractions, so segment colours (Okabe-Ito) are distinguishable.
LINEAGE = {
    'B cells':        ['B cells naive', 'B cells memory'],
    'Plasma cells':   ['Plasma cells'],
    'CD8 T':          ['T cells CD8'],
    'CD4 T':          ['T cells CD4 naive', 'T cells CD4 memory resting',
                       'T cells CD4 memory activated', 'T cells follicular helper',
                       'T cells regulatory (Tregs)'],
    'gamma-delta T':  ['T cells gamma delta'],
    'NK':             ['NK cells resting', 'NK cells activated'],
    'Monocytes':      ['Monocytes'],
    'Macrophages':    ['Macrophages M0', 'Macrophages M1', 'Macrophages M2'],
    'Dendritic':      ['Dendritic cells resting', 'Dendritic cells activated'],
    'Mast':           ['Mast cells resting', 'Mast cells activated'],
    'Eosinophils':    ['Eosinophils'],
    'Neutrophils':    ['Neutrophils'],
}
lineage = pd.DataFrame(
    {name: frac[[c for c in cols if c in frac.columns]].sum(axis=1)
     for name, cols in LINEAGE.items()},
    index=frac.index)

order = labels.sort_values().index               # CON | CSR | DCM contiguous
fr = lineage.loc[order]
g = labels.loc[order].values
fig, ax = plt.subplots(figsize=(7.0, 3.4))
# 12 distinguishable hues (Okabe-Ito core + Tol-derived extras) so no two lineages
# collide; F.OKABE alone (9) would recycle green onto both CD8 T and Neutrophils.
# Keyed by lineage name so the dominant ~50% Neutrophil segment gets a soft neutral
# (not heavy black), and dark colours go to thin segments.
LIN_COLORMAP = {
    'B cells':       '#0072B2',  # blue
    'Plasma cells':  '#D55E00',  # vermillion
    'CD8 T':         '#009E73',  # bluish-green
    'CD4 T':         '#CC79A7',  # reddish-purple
    'gamma-delta T': '#E69F00',  # orange
    'NK':            '#56B4E9',  # sky blue
    'Monocytes':     '#F0E442',  # yellow
    'Macrophages':   '#332288',  # indigo
    'Dendritic':     '#882255',  # wine
    'Mast':          '#44AA99',  # teal
    'Eosinophils':   '#999933',  # olive
    'Neutrophils':   '#D9D9D9',  # light grey (dominant segment -> neutral, not black)
}
cols_cycle = [LIN_COLORMAP[c] for c in fr.columns]
bottom = np.zeros(len(order))
x = np.arange(len(order))
for i, c in enumerate(fr.columns):
    ax.bar(x, fr[c].values, bottom=bottom, color=cols_cycle[i], width=1.0,
           label=c, edgecolor='white', linewidth=0.3)
    bottom += fr[c].values
# labelled group spans (replaces the 2 thin black separator lines)
ax.set_xlim(-0.5, len(order)-0.5); ax.set_ylim(0, 1.12)
span_shade = {'CON': '#F2F2F2', 'CSR': '#FFFFFF', 'DCM': '#F2F2F2'}
counts = {grp: int((g == grp).sum()) for grp in ['CON', 'CSR', 'DCM']}
start = 0
for grp in ['CON', 'CSR', 'DCM']:
    idx = np.where(g == grp)[0]
    if len(idx) == 0:
        continue
    lo, hi = idx.min() - 0.5, idx.max() + 0.5
    ax.axvspan(lo, hi, ymin=0, ymax=1.0/1.12, color=span_shade[grp], zorder=0)
    ax.text((lo + hi) / 2, 1.04, f'{grp} (n={counts[grp]})',
            ha='center', va='bottom', fontsize=8, fontweight='bold')
ax.set_ylabel('Estimated immune-cell fraction')
ax.set_xlabel('Samples (CON | CSR | DCM)')
ax.set_xticks([])
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
sns.despine(ax=ax)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=5,
          fontsize=7, frameon=False, handlelength=1.0, columnspacing=1.2,
          handletextpad=0.5)
ax.set_title('Blood immune-cell composition (CIBERSORT)')
save(fig, 'Fig_immune_stackedbar')

# ---- Fig 2: boxplots of immune cells differing across groups ----
def mwu(a, b):
    a, b = a.dropna(), b.dropna()
    if len(a) < 2 or len(b) < 2 or (a.std()==0 and b.std()==0): return np.nan
    try: return stats.mannwhitneyu(a, b).pvalue
    except Exception: return np.nan

df = frac.copy(); df['group'] = labels
stat_rows = []
for c in frac.columns:
    p_cerv = mwu(df.loc[labels.isin(['CSR','DCM']), c], df.loc[labels=='CON', c])
    p_csr  = mwu(df.loc[labels=='CSR', c], df.loc[labels=='CON', c])
    stat_rows.append([c, frac[c].mean(), p_cerv, p_csr])
cell_stats = pd.DataFrame(stat_rows, columns=['cell','mean_frac','p_CERVvsCON','p_CSRvsCON']).sort_values('p_CERVvsCON')
cell_stats.to_csv(os.path.join(ROOT, 'immune_group_stats.csv'), index=False)
sig_cells = cell_stats[(cell_stats.p_CERVvsCON < 0.05) | (cell_stats.p_CSRvsCON < 0.05)]['cell'].tolist()
plot_cells = sig_cells if len(sig_cells) >= 4 else cell_stats.head(8)['cell'].tolist()
n = len(plot_cells); ncol = 4; nrow = int(np.ceil(n/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(min(7.0, 1.75*ncol), 2.5*nrow)); axes = np.array(axes).ravel()
for k, (ax, c) in enumerate(zip(axes, plot_cells)):
    F.boxstrip(ax, df, 'group', c, ['CON','CSR','DCM'], PAL)
    # p-value pools CSR+DCM vs CON; comparison stated once in the figure note below
    p = cell_stats.set_index('cell').loc[c, 'p_CERVvsCON']
    s = F.star(p)
    ax.set_title(f'{c.replace("_"," ")}\nP = {p:.3f}{("  "+s) if s else ""}', fontsize=7.5)
    ax.set_ylabel('Fraction')
    sns.despine(ax=ax)
    F.panel(ax, chr(ord('A')+k))
for ax in axes[n:]: ax.axis('off')
# state the tested comparison once so P-values match the displayed groups
fig.text(0.5, -0.02, 'P: Cervical degeneration (CSR + DCM, n=39) vs CON (n=7), Mann-Whitney U',
         ha='center', va='top', fontsize=7)
fig.tight_layout(); save(fig, 'Fig_immune_boxplots')

# ================= 2. Hub-gene ROC (blood-expressed) =================
hub_blood = [h for h in ['IL1B','MMP9','AKT1','TP53'] if h in log2.index]
def roc_ci(y, score, n_boot=2000, seed=0):
    rng = np.random.RandomState(seed)
    fpr, tpr, _ = roc_curve(y, score); a = auc(fpr, tpr)
    if a < 0.5: fpr, tpr, _ = roc_curve(y, -score); a = auc(fpr, tpr); score = -score
    aucs = []
    idx = np.arange(len(y))
    for _ in range(n_boot):
        bi = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(y[bi])) < 2: continue
        f2, t2, _ = roc_curve(y[bi], score[bi]); aucs.append(auc(f2, t2))
    lo, hi = np.percentile(aucs, [2.5, 97.5]) if aucs else (np.nan, np.nan)
    return fpr, tpr, a, lo, hi

CURVE_LS = ['-', '--', '-.', ':', '-']   # linestyle redundancy so curves separate in grayscale/where they cross
ROC_TITLE = {'CERVvsCON': 'ROC: Cerv. degen. vs CON',
             'CSRvsCON':  'ROC: CSR vs CON'}
for contrast, mask_case in [('CERVvsCON', labels.isin(['CSR','DCM'])), ('CSRvsCON', labels=='CSR')]:
    sel = mask_case | (labels=='CON')
    y = (labels[sel]=='CON').map({True:0}).reindex(labels[sel].index)
    y = mask_case[sel].astype(int).values  # case=1, con=0
    X = log2.loc[hub_blood, sel.index[sel.values]].T
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    rows=[]
    for j, h in enumerate(hub_blood):
        fpr,tpr,a,lo,hi = roc_ci(y, X[h].values)
        ax.plot(fpr, tpr, lw=1.8, color=F.OKABE[j % len(F.OKABE)],
                ls=CURVE_LS[j % len(CURVE_LS)],
                label=f'{h}: AUC={a:.2f} ({lo:.2f}-{hi:.2f})'); rows.append([h,a,lo,hi])
    # combined logistic (LOOCV)
    if X.shape[0] >= 8:
        proba = np.zeros(len(y))
        for tr, te in LeaveOneOut().split(X):
            m = LogisticRegression(max_iter=1000).fit(X.iloc[tr], y[tr]); proba[te]=m.predict_proba(X.iloc[te])[:,1]
        fpr,tpr,a,lo,hi = roc_ci(y, proba)
        ax.plot(fpr, tpr, color='k', ls=(0, (3, 1, 1, 1)), lw=2.0,
                label=f'Combined (LOOCV): AUC={a:.2f} ({lo:.2f}-{hi:.2f})'); rows.append(['Combined_LOOCV',a,lo,hi])
    ax.plot([0,1],[0,1],color='grey',ls=':',lw=0.9)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02); ax.set_aspect('equal')
    ax.set_xlabel('1 - Specificity'); ax.set_ylabel('Sensitivity')
    ax.set_title(ROC_TITLE[contrast])
    sns.despine(ax=ax)
    ax.legend(fontsize=6.5, loc='lower right', frameon=False)
    save(fig, f'Fig_ROC_{contrast}')
    pd.DataFrame(rows,columns=['feature','AUC','CI_low','CI_high']).to_csv(os.path.join(ROOT,f'ROC_{contrast}.csv'),index=False)
    print(f'ROC {contrast}: n_case={int(y.sum())} n_con={int((y==0).sum())}')

# ---- Fig 3: hub-gene expression boxplots (2x2 grid, shared y-label) ----
nh = len(hub_blood); ncol = 2; nrow = int(np.ceil(nh/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.6, 4.4)); axes = np.array(axes).ravel()
hd = log2.loc[hub_blood].T.copy(); hd['group'] = labels
for k, (ax, h) in enumerate(zip(axes, hub_blood)):
    F.boxstrip(ax, hd, 'group', h, ['CON','CSR','DCM'], PAL)
    # p-value pools CSR+DCM (cervical degeneration) vs CON; comparison stated once below
    p = mwu(hd.loc[labels.isin(['CSR','DCM']), h], hd.loc[labels=='CON', h])
    s = F.star(p)
    ax.set_title(f'{h}\nP = {p:.3f}{("  "+s) if s else ""}', fontsize=7.5)
    ax.set_ylabel('')
    sns.despine(ax=ax)
    F.panel(ax, chr(ord('A')+k))
for ax in axes[nh:]: ax.axis('off')
fig.supylabel('log2(FPKM+1)', fontsize=8)
# state the tested comparison once so P-values match the displayed groups
fig.text(0.5, -0.03, 'P: Cervical degeneration (CSR + DCM, n=39) vs CON (n=7), Mann-Whitney U',
         ha='center', va='top', fontsize=7)
fig.tight_layout(); save(fig, 'Fig_hubgene_expr')

# ================= 3. hub gene x immune correlation =================
corr = pd.DataFrame(index=hub_blood, columns=frac.columns, dtype=float)
pmat = corr.copy()
for h in hub_blood:
    for c in frac.columns:
        r,p = stats.spearmanr(log2.loc[h], frac[c])
        corr.loc[h,c]=r; pmat.loc[h,c]=p
corr.to_csv(os.path.join(ROOT,'hub_immune_spearman_r.csv'))
# Transpose: 22 immune cells as horizontal y-labels, 4 hub genes on x. Contiguous matrix
# (no blank spacer columns); de-underscore labels; square cells; stroked-white black asterisks.
corr_t = corr.astype(float).T            # immune cells (rows) x hub genes (cols)
pmat_t = pmat.astype(float).T
corr_t.index = [str(c).replace('_', ' ') for c in corr_t.index]
fig, ax = plt.subplots(figsize=(3.8, 6.4))
sns.heatmap(corr_t, cmap=F.DIVERGING, center=0, vmin=-.7, vmax=.7,
            annot=False, cbar_kws={'label': 'Spearman r', 'shrink': 0.5},
            ax=ax, linewidths=0.4, linecolor='white', square=True)
# black asterisks with white stroke read on any cell colour
for i in range(corr_t.shape[0]):
    for j in range(corr_t.shape[1]):
        if pmat_t.iloc[i, j] < 0.05:
            t = ax.text(j + 0.5, i + 0.5, '*', ha='center', va='center',
                        color='k', fontsize=9, fontweight='bold')
            F.stroke(t)
ax.set_xlabel('Hub gene'); ax.set_ylabel('Immune cell')
ax.set_title('Hub gene-immune correlation')
plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7.5)
plt.setp(ax.get_yticklabels(), rotation=0, fontsize=7)
fig.tight_layout(); save(fig, 'Fig_hub_immune_corr')

print('\n==== SUMMARY ====')
print('Top immune cells (CERV vs CON):'); print(cell_stats.head(6).to_string(index=False))
print('\nFigures + tables written to', FIG)
