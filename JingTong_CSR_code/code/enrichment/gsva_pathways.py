# -*- coding: utf-8 -*-
"""
GSVA pathway-activity scoring on GSE223227: score every sample for KEGG pathways, then
compare the network-pharmacology pathways between cervical degeneration (CSR+DCM) and control.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
from scipy import stats
import gseapy as gp

ML=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'ML_GSE223227')
OUT=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(OUT,'figures')
expr=pd.read_csv(os.path.join(ML,'expr_symbol.csv.gz'),index_col=0)
labels=pd.read_csv(os.path.join(ML,'labels.csv'),index_col=0)['group']
log2=np.log2(expr+1); log2=log2[log2.mean(1)>0.5]
PAL={'CON':'#4C72B0','CSR':'#C44E52','DCM':'#DD8452'}

gv=gp.gsva(data=log2, gene_sets='KEGG_2021_Human', min_size=10, max_size=500, threads=4, outdir=None)
es=gv.res2d.copy()
# pivot to pathway x sample
val=[c for c in es.columns if c.lower() in ('es','nes','score')]
mat=es.pivot(index='Term',columns='Name',values=val[0]).astype(float)[labels.index]
mat.to_csv(os.path.join(OUT,'GSVA_kegg_scores.csv'))

KEY=['PI3K-Akt','TNF signaling','IL-17','MAPK signaling','AGE-RAGE','NF-kappa B',
     'Osteoclast','Chemokine','Lipid and atherosclerosis','Cytokine-cytokine']
def mwu(a,b):
    try:return stats.mannwhitneyu(a.dropna(),b.dropna()).pvalue
    except Exception:return np.nan
rows=[]
sel=[]
for t in mat.index:
    if any(k.lower() in t.lower() for k in KEY):
        p=mwu(mat.loc[t,labels.isin(['CSR','DCM'])], mat.loc[t,labels=='CON'])
        d=mat.loc[t,labels.isin(['CSR','DCM'])].mean()-mat.loc[t,labels=='CON'].mean()
        rows.append([t,round(d,3),round(p,4)]); sel.append(t)
res=pd.DataFrame(rows,columns=['pathway','delta_CERVvsCON','p']).sort_values('p')
res.to_csv(os.path.join(OUT,'GSVA_key_pathways.csv'),index=False)
print(res.to_string(index=False))

if sel:
    long=mat.loc[sel].T.copy(); long['group']=labels
    m=long.melt(id_vars='group',var_name='pathway',value_name='GSVA')
    m['short']=m['pathway'].str.replace(r' \(.*','',regex=True).str.slice(0,28)
    fig,ax=plt.subplots(figsize=(min(14,1.3*len(sel)+3),5))
    sns.boxplot(data=m,x='short',y='GSVA',hue='group',hue_order=['CON','CSR','DCM'],palette=PAL,ax=ax,fliersize=0)
    plt.setp(ax.get_xticklabels(),rotation=35,ha='right',fontsize=8)
    ax.set_title('GSVA pathway activity by group (key network-pharmacology pathways)')
    ax.set_xlabel(''); ax.legend(fontsize=8,frameon=False); fig.tight_layout()
    for ext in ('png','pdf'): fig.savefig(os.path.join(FIG,'Fig_GSVA_pathways.'+ext),dpi=300,bbox_inches='tight')
    print('\nSaved Fig_GSVA_pathways + GSVA_*.csv')
