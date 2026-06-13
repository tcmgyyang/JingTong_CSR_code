# -*- coding: utf-8 -*-
"""
Task 1 robustness for the blood module (GSE223227):
 (A) Age-confound transparency: age distributions + hub-gene~age correlation
     (CON 23-28 vs patients 33-69 do NOT overlap -> age adjustment is not identifiable;
      we therefore QUANTIFY and DISCLOSE the confound honestly).
 (B) ssGSEA immune cross-check with a compact canonical immune-marker panel, compared
     against the CIBERSORT result, plus correlation with blood hub genes IL1B/MMP9.
"""
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F
F.apply()

ROOT=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(ROOT,'figures')
expr=pd.read_csv(os.path.join(ROOT,'expr_symbol.csv.gz'),index_col=0)
labels=pd.read_csv(os.path.join(ROOT,'labels.csv'),index_col=0)['group']
meta=pd.read_csv(os.path.join(ROOT,'sample_meta.csv'),index_col=0)
expr=expr[labels.index]; log2=np.log2(expr+1)
meta=meta.loc[labels.index]
PAL=F.PAL  # Okabe-Ito colourblind-safe triad: CON #0072B2, CSR #D55E00, DCM #009E73

# ---------- (A) age confound ----------
hub_blood=['IL1B','MMP9','AKT1','TP53']
agecorr=[]
for h in hub_blood:
    r,p=stats.spearmanr(meta['age'],log2.loc[h]); agecorr.append([h,round(r,2),round(p,4)])
tbl=pd.DataFrame(agecorr,columns=['hub_gene','Spearman r vs age','p'])
tbl.to_csv(os.path.join(ROOT,'hubgene_age_correlation.csv'),index=False)

fig,axes=plt.subplots(1,2,figsize=(7.0,3.2),gridspec_kw={'width_ratios':[1,1.05]})
# A: age distribution by group (Okabe-Ito boxes + jittered points)
dfa=meta.assign(group=labels)
sns.boxplot(data=dfa,x='group',y='age',order=['CON','CSR','DCM'],hue='group',
            hue_order=['CON','CSR','DCM'],palette=PAL,legend=False,
            ax=axes[0],fliersize=0,linewidth=0.8,width=0.62)
sns.stripplot(data=dfa,x='group',y='age',order=['CON','CSR','DCM'],color='0.15',
              size=3,alpha=0.8,jitter=0.18,ax=axes[0])
axes[0].set_title('Age by group'); axes[0].set_ylabel('Age (years)'); axes[0].set_xlabel('')
F.panel(axes[0],'A')

# B: horizontal lollipop of hub-gene ~ age Spearman r (stem + dot), values to 2 dp
ax=axes[1]
order=tbl.iloc[::-1].reset_index(drop=True)        # top gene at top of axis
ypos=np.arange(len(order))
rvals=order['Spearman r vs age'].astype(float).values
ax.axvline(0,color=F.NEUTRAL,lw=0.8,zorder=0)
ax.hlines(ypos,0,rvals,color='0.55',lw=1.4,zorder=1)
ax.scatter(rvals,ypos,s=55,color=F.OKABE[0],edgecolor='white',linewidth=0.8,zorder=2)
for y,r in zip(ypos,rvals):
    off=0.015 if r>=0 else -0.015
    ax.text(r+off,y,f'{r:.2f}',va='center',ha='left' if r>=0 else 'right',fontsize=7)
ax.set_yticks(ypos); ax.set_yticklabels(order['hub_gene'])
ax.set_xlabel('Spearman r vs age'); ax.set_title('Hub-gene correlation with age')
xpad=max(0.1,np.abs(rvals).max()*1.35)
ax.set_xlim(-xpad,xpad); ax.set_ylim(-0.5,len(order)-0.5)
sns.despine(ax=ax,left=True); ax.tick_params(left=False)
F.panel(ax,'B')
fig.tight_layout()
F.save(fig,os.path.join(FIG,'Fig_age_confound'))
print('Age confound:'); print(meta.assign(group=labels).groupby('group')['age'].agg(['mean','min','max']).round(1).to_string())
print('\nHub-gene vs age:'); print(tbl.to_string(index=False))

# ---------- (B) ssGSEA cross-check ----------
IMMUNE={
 'CD8_T_cell':['CD8A','CD8B','GZMK','GZMA','NKG7','CD3D'],
 'CD4_T_cell':['CD4','IL7R','CCR7','CD3D','CD3E','TCF7'],
 'Treg':['FOXP3','IL2RA','CTLA4','IKZF2'],
 'B_cell':['CD19','MS4A1','CD79A','CD79B','TNFRSF13B'],
 'Plasma_cell':['MZB1','IGHG1','JCHAIN','XBP1','CD38'],
 'NK_cell':['KLRD1','NCR1','GNLY','KLRF1','NCAM1'],
 'Monocyte':['CD14','LYZ','FCN1','S100A8','S100A9'],
 'Macrophage':['CD68','CD163','MRC1','CSF1R'],
 'Dendritic_cell':['ITGAX','CD1C','CLEC9A','LILRA4'],
 'Neutrophil':['FCGR3B','CSF3R','CEACAM3','FUT4','S100A12'],
 'Mast_cell':['TPSAB1','TPSB2','CPA3','MS4A2'],
 'Cytotoxic':['GZMB','PRF1','GNLY','GZMA'],
}
ss=None
try:
    import gseapy as gp
    r=gp.ssgsea(data=log2, gene_sets=IMMUNE, sample_norm_method='rank', outdir=None, threads=4, min_size=3)
    d=r.res2d.copy()
    val='NES' if 'NES' in d.columns else ('ES' if 'ES' in d.columns else d.columns[-1])
    ss=d.pivot(index='Name',columns='Term',values=val).astype(float)
    print('\nssGSEA via gseapy OK', ss.shape)
except Exception as e:
    print('\ngseapy ssGSEA failed (', e, ') -> fallback: mean z-scored marker score')
    z=log2.sub(log2.mean(1),axis=0).div(log2.std(1)+1e-9,axis=0)
    ss=pd.DataFrame({c:z.reindex([g for g in genes if g in z.index]).mean() for c,genes in IMMUNE.items()})
ss=ss.reindex(labels.index)
ss.to_csv(os.path.join(ROOT,'ssGSEA_immune.csv'))

# de-underscore cell labels for display (CD8_T_cell -> CD8 T cell)
def nice(c): return c.replace('_',' ')

# group comparison: combined cervical patients (CSR+DCM) vs CON
def mwu(a,b):
    a,b=a.dropna(),b.dropna()
    try:return stats.mannwhitneyu(a,b).pvalue
    except Exception:return np.nan
rows=[]; pmap={}
for c in ss.columns:
    p=mwu(ss.loc[labels.isin(['CSR','DCM']),c],ss.loc[labels=='CON',c])
    rows.append([c,round(p,3)]); pmap[c]=p
gp_stats=pd.DataFrame(rows,columns=['cell','p_CERVvsCON']).sort_values('p_CERVvsCON')
gp_stats.to_csv(os.path.join(ROOT,'ssGSEA_group_stats.csv'),index=False)

fig,ax=plt.subplots(figsize=(7.0,3.6))
cells=list(ss.columns)
dd=ss.copy(); dd['group']=labels
m=dd.melt(id_vars='group',var_name='cell',value_name='score')
m['cell']=m['cell'].map(nice)
order_nice=[nice(c) for c in cells]
sns.boxplot(data=m,x='cell',y='score',order=order_nice,hue='group',hue_order=['CON','CSR','DCM'],
            palette=PAL,ax=ax,fliersize=0,linewidth=0.6,width=0.74)
# significance brackets: cervical (CSR+DCM) vs CON, F.star convention
ymax=m['score'].max(); yr=m['score'].max()-m['score'].min()
ytop=ymax+0.04*yr
for i,c in enumerate(cells):
    s=F.star(pmap[c])
    if not s: continue
    t=ax.text(i,ytop,s,ha='center',va='bottom',fontsize=8,fontweight='bold',color='0.15')
ax.set_ylim(top=ytop+0.12*yr)
plt.setp(ax.get_xticklabels(),rotation=40,ha='right')
ax.set_xlabel(''); ax.set_ylabel('ssGSEA enrichment score')
ax.set_title('ssGSEA immune scores by group')
ax.legend(title='',loc='upper right',frameon=False,ncol=3,handlelength=1.1,columnspacing=1.0)
ax.text(0.0,1.04,'* cervical (CSR+DCM) vs CON, Mann-Whitney',transform=ax.transAxes,
        fontsize=6.5,color='0.4',ha='left',va='bottom')
fig.tight_layout()
F.save(fig,os.path.join(FIG,'Fig_ssGSEA_immune'))

# correlation hub genes x ssGSEA cells
corr=pd.DataFrame(index=hub_blood,columns=ss.columns,dtype=float); pm=corr.copy()
for h in hub_blood:
    for c in ss.columns:
        r,p=stats.spearmanr(log2.loc[h],ss[c]); corr.loc[h,c]=r; pm.loc[h,c]=p
cval=corr.astype(float)
vmax=float(np.nanmax(np.abs(cval.values)))      # symmetric scale, center=0
fig,ax=plt.subplots(figsize=(7.0,3.0))
sns.heatmap(cval,cmap=F.DIVERGING,center=0,vmin=-vmax,vmax=vmax,square=True,
            ax=ax,linewidths=.5,linecolor='white',
            xticklabels=[nice(c) for c in cval.columns],yticklabels=list(cval.index),
            cbar_kws={'label':'Spearman r','shrink':0.75})
# asterisks where p<0.05, black glyph with white stroke for legibility on any cell
for i,h in enumerate(corr.index):
    for j,c in enumerate(corr.columns):
        if pm.loc[h,c]<0.05:
            t=ax.text(j+0.5,i+0.5,'*',ha='center',va='center',fontsize=10,
                      fontweight='bold',color='black')
            F.stroke(t)
ax.set_xlabel('Immune cell type'); ax.set_ylabel('Hub gene')
ax.set_title('Hub genes vs ssGSEA immune scores')
plt.setp(ax.get_xticklabels(),rotation=40,ha='right')
plt.setp(ax.get_yticklabels(),rotation=0)
fig.tight_layout()
F.save(fig,os.path.join(FIG,'Fig_ssGSEA_hubcorr'))
print('\nssGSEA top group diffs:'); print(gp_stats.head(6).to_string(index=False))
print('\nFigures: Fig_age_confound, Fig_ssGSEA_immune, Fig_ssGSEA_hubcorr')
