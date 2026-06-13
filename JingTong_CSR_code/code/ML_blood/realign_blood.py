# -*- coding: utf-8 -*-
"""
Realign the blood diagnostic module to the 4-source multi-database hub set.
Blood-expressed hubs (FPKM>=1 in GSE223227) among the 18 >=2-source hubs:
  IL1B, NFKBIA, HIF1A, BCL2L1, CXCL8, MMP9, TNF, AKT1, MYC, TP53.
Produces: per-gene diagnostic AUC (+95% CI) lollipop; LASSO-selected combined panel ROC
(LOOCV); hub-gene expression boxplots; hub gene-immune correlation.
"""
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F; F.apply()
ROOT = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(ROOT, 'figures')

expr = pd.read_csv(os.path.join(ROOT, 'expr_symbol.csv.gz'), index_col=0)
labels = pd.read_csv(os.path.join(ROOT, 'labels.csv'), index_col=0)['group']
expr = expr[labels.index]; log2 = np.log2(expr + 1)
PAL = F.PAL
HUB18 = ['AKT1','IL6','TNF','IL1B','ESR1','MYC','EGFR','JUN','HIF1A','MMP9','CCL2',
         'CXCL8','BDNF','IL10','IL2','IL4','BCL2L1','NFKBIA']
hub_blood = [h for h in HUB18 if h in log2.index and expr.loc[h].mean() >= 1]
print('blood-expressed hubs:', hub_blood)

def roc_ci(y, score, nb=2000, seed=0):
    rng = np.random.RandomState(seed)
    fpr, tpr, _ = roc_curve(y, score); a = auc(fpr, tpr)
    if a < 0.5: score = -score; fpr, tpr, _ = roc_curve(y, score); a = auc(fpr, tpr)
    aucs = []
    for _ in range(nb):
        bi = rng.choice(len(y), len(y), replace=True)
        if len(np.unique(y[bi])) < 2: continue
        f, t, _ = roc_curve(y[bi], score[bi]); aucs.append(auc(f, t))
    lo, hi = np.percentile(aucs, [2.5, 97.5]) if aucs else (np.nan, np.nan)
    return fpr, tpr, a, lo, hi

CONTRASTS = [('Cervical degeneration vs CON', labels.isin(['CSR','DCM'])),
             ('CSR vs CON', labels == 'CSR')]
auc_tab = {}
combined = {}
for cname, mask in CONTRASTS:
    sel = (mask | (labels == 'CON'))
    y = mask[sel].astype(int).values
    X = log2.loc[hub_blood, sel.index[sel.values]].T
    rows = []
    for h in hub_blood:
        _, _, a, lo, hi = roc_ci(y, X[h].values)
        rows.append([h, a, lo, hi])
    auc_tab[cname] = pd.DataFrame(rows, columns=['gene','AUC','lo','hi']).sort_values('AUC')
    # LASSO-selected combined panel (L1 logistic, CV lambda) then LOOCV AUC
    Xs = StandardScaler().fit_transform(X.values)
    try:
        las = LogisticRegressionCV(Cs=20, cv=min(5, y.sum(), (y==0).sum()), penalty='l1',
                                   solver='liblinear', scoring='roc_auc', max_iter=2000).fit(Xs, y)
        keep = [hub_blood[i] for i in np.where(np.abs(las.coef_[0]) > 1e-6)[0]]
    except Exception:
        keep = hub_blood
    keep = keep or hub_blood
    Xk = X[keep].values
    proba = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(Xk):
        sc = StandardScaler().fit(Xk[tr])
        m = LogisticRegression(max_iter=2000).fit(sc.transform(Xk[tr]), y[tr])
        proba[te] = m.predict_proba(sc.transform(Xk[te]))[:, 1]
    fpr, tpr, a, lo, hi = roc_ci(y, proba)
    combined[cname] = (fpr, tpr, a, lo, hi, keep)
    print(f'{cname}: LASSO panel={keep} combined AUC={a:.2f} ({lo:.2f}-{hi:.2f})')

pd.concat({k: v.set_index('gene') for k, v in auc_tab.items()}, axis=1).to_csv(os.path.join(ROOT,'blood_hub_AUC.csv'))

# ---- Fig A: per-gene AUC lollipop (both contrasts) + combined ROC ----
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.2), gridspec_kw={'width_ratios':[1,1,1.1]})
for ax, (cname, _) in zip(axes[:2], CONTRASTS):
    t = auc_tab[cname]; ypos = np.arange(len(t))
    sig = t['lo'] > 0.5
    ax.hlines(ypos, 0.5, t['AUC'], color='#CCCCCC', lw=1.2, zorder=1)
    ax.errorbar(t['AUC'], ypos, xerr=[t['AUC']-t['lo'], t['hi']-t['AUC']], fmt='none',
                ecolor='#999999', elinewidth=0.8, capsize=2, zorder=2)
    ax.scatter(t['AUC'], ypos, s=34, c=[F.PAL['CSR'] if s else F.PAL['CON'] for s in sig], zorder=3, edgecolors='black', linewidths=0.4)
    ax.set_yticks(ypos); ax.set_yticklabels(t['gene'], fontsize=7)
    ax.axvline(0.5, color='k', ls=':', lw=0.8); ax.set_xlim(0.3, 1.02)
    ax.set_xlabel('AUC (95% CI)'); ax.set_title(cname, fontsize=8)
axA = axes[2]
for cname, ls in zip([c[0] for c in CONTRASTS], ['-', '--']):
    fpr, tpr, a, lo, hi, keep = combined[cname]
    axA.plot(fpr, tpr, ls=ls, lw=1.8, color=F.PAL['DCM'],
             label=f'{cname.split(" vs")[0][:10]}\nAUC={a:.2f}')
axA.plot([0,1],[0,1],color='grey',ls=':',lw=0.8); axA.set_aspect('equal')
axA.set_xlabel('1 - Specificity'); axA.set_ylabel('Sensitivity')
axA.set_title('LASSO panel (LOOCV)', fontsize=8); axA.legend(fontsize=6.5, loc='lower right')
for i, ax in enumerate(axes): F.panel(ax, chr(65+i))
fig.suptitle('Blood diagnostic value of multi-source hub genes (GSE223227)', fontsize=10)
fig.tight_layout(); F.save(fig, os.path.join(FIG, 'Fig_blood_hub_AUC'))

# ---- Fig B: expression boxplots (2x5) ----
nr, nc = 2, 5
fig, axes = plt.subplots(nr, nc, figsize=(7.2, 3.6)); axes = np.array(axes).ravel()
hd = log2.loc[hub_blood].T.copy(); hd['group'] = labels
def mwu(a, b):
    try: return stats.mannwhitneyu(a.dropna(), b.dropna()).pvalue
    except Exception: return np.nan
for k, (ax, h) in enumerate(zip(axes, hub_blood)):
    F.boxstrip(ax, hd, 'group', h, ['CON','CSR','DCM'], PAL)
    p = mwu(hd.loc[labels.isin(['CSR','DCM']), h], hd.loc[labels=='CON', h]); s = F.star(p)
    ax.set_title(f'{h}\nP={p:.3f}{("  "+s) if s else ""}', fontsize=7); ax.set_ylabel('')
    F.panel(ax, chr(65+k))
for ax in axes[len(hub_blood):]: ax.axis('off')
fig.supylabel('log2(FPKM+1)', fontsize=8)
fig.text(0.5, -0.02, 'P: cervical degeneration (CSR+DCM, n=39) vs CON (n=7), Mann-Whitney U', ha='center', va='top', fontsize=7)
fig.tight_layout(); F.save(fig, os.path.join(FIG, 'Fig_hubgene_expr'))

# ---- Fig C: hub gene x immune correlation (CIBERSORT) ----
frac = pd.read_csv(os.path.join(ROOT, 'CIBERSORT_fractions.csv'), index_col=0).reindex(labels.index)
corr = pd.DataFrame(index=hub_blood, columns=frac.columns, dtype=float); pm = corr.copy()
for h in hub_blood:
    for c in frac.columns:
        r, p = stats.spearmanr(log2.loc[h], frac[c]); corr.loc[h,c]=r; pm.loc[h,c]=p
ct = corr.astype(float).T; ct.index=[str(c).replace('_',' ') for c in ct.index]
vmax = float(np.nanmax(np.abs(ct.values)))
fig, ax = plt.subplots(figsize=(5.4, 6.6))
sns.heatmap(ct, cmap=F.DIVERGING, center=0, vmin=-vmax, vmax=vmax, ax=ax, square=True,
            linewidths=0.4, linecolor='white', cbar_kws={'label':'Spearman r','shrink':0.5})
for i, cell in enumerate(ct.index):
    for j, g in enumerate(ct.columns):
        if pm.T.values[i, j] < 0.05:
            t = ax.text(j+0.5, i+0.5, '*', ha='center', va='center', fontsize=9, color='black'); F.stroke(t)
ax.set_xlabel('Hub gene'); ax.set_ylabel('Immune cell'); plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
ax.set_title('Hub gene - immune-cell correlation (* p<0.05)', fontsize=9)
fig.tight_layout(); F.save(fig, os.path.join(FIG, 'Fig_hub_immune_corr'))
print('\nSaved Fig_blood_hub_AUC, Fig_hubgene_expr (10 genes), Fig_hub_immune_corr (10 genes)')
