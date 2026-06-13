# -*- coding: utf-8 -*-
"""
Single-gene GSEA: for each diagnostic hub gene, rank all genes by their correlation with
that gene across GSE223227, then GSEA (KEGG) -> the pathways each biomarker drives.
This links the ML-selected biomarkers back to the network-pharmacology pathways.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp

ML=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'ML_GSE223227')
FIG=os.path.join(os.path.dirname(os.path.abspath(__file__)),'figures'); os.makedirs(FIG,exist_ok=True)
expr=pd.read_csv(os.path.join(ML,'expr_symbol.csv.gz'),index_col=0)
log2=np.log2(expr+1)
# drop low-expression genes (noise) for cleaner ranking
log2=log2[log2.mean(1)>0.5]
ranks=log2.rank(axis=1)                      # for Spearman = Pearson on ranks
GENES=['AKT1','TNF','IL1B','MYC','HIF1A','MMP9','CXCL8','BCL2L1','NFKBIA']  # 9 blood-expressed multi-source hubs
KEY_KEGG=['PI3K-Akt','TNF','IL-17','MAPK','AGE-RAGE','NF-kappa B','Lipid and atherosclerosis','TNF signaling']

summary=[]
for g in GENES:
    if g not in ranks.index: print('skip',g); continue
    gr=ranks.loc[g]
    # Pearson corr of every gene's ranks with target gene's ranks = Spearman
    R=ranks.sub(ranks.mean(1),axis=0)
    num=R.mul(gr-gr.mean(),axis=1).sum(1)
    den=np.sqrt((R**2).sum(1))*np.sqrt(((gr-gr.mean())**2).sum())
    rho=(num/den).drop(g).dropna().sort_values(ascending=False)
    rnk=rho.reset_index(); rnk.columns=['gene','score']
    try:
        pre=gp.prerank(rnk=rnk, gene_sets='KEGG_2021_Human', min_size=10, max_size=500,
                       permutation_num=1000, seed=42, threads=4, outdir=None, no_plot=True)
        res=pre.res2d.copy()
        res['NES']=pd.to_numeric(res['NES'],errors='coerce')
        res['FDR q-val']=pd.to_numeric(res['FDR q-val'],errors='coerce')
        res=res.sort_values('NES',ascending=False)
        res.to_csv(os.path.join(os.path.dirname(FIG),f'GSEA_{g}.csv'),index=False)
        top=res[res['FDR q-val']<0.25].head(8)
        print(f'\n=== {g}: top enriched KEGG (FDR<0.25) ===')
        print(top[['Term','NES','FDR q-val']].to_string(index=False))
        for _,row in res.iterrows():
            if any(k.lower() in str(row['Term']).lower() for k in KEY_KEGG):
                summary.append([g,row['Term'],row['NES'],row['FDR q-val']])
    except Exception as e:
        print(f'{g}: GSEA failed -> {e}')

# heatmap: NES of key network-pharmacology pathways across the biomarkers
if summary:
    s=pd.DataFrame(summary,columns=['gene','Term','NES','FDR'])
    piv=s.pivot_table(index='Term',columns='gene',values='NES',aggfunc='first')
    import seaborn as sns
    fig,ax=plt.subplots(figsize=(7,max(3,0.5*len(piv))))
    sns.heatmap(piv.astype(float),cmap='RdBu_r',center=0,annot=True,fmt='.2f',ax=ax,
                cbar_kws={'label':'GSEA NES'},linewidths=.4)
    ax.set_title('Single-gene GSEA: network-pharmacology pathways driven by each biomarker')
    plt.setp(ax.get_yticklabels(),fontsize=8); fig.tight_layout()
    for ext in ('png','pdf'): fig.savefig(os.path.join(FIG,'Fig_single_gene_GSEA.'+ext),dpi=300,bbox_inches='tight')
    s.to_csv(os.path.join(os.path.dirname(FIG),'GSEA_key_pathways_summary.csv'),index=False)
    print('\nSaved Fig_single_gene_GSEA + GSEA_*.csv')
else:
    print('\nNo key-pathway hits to plot (see per-gene GSEA_*.csv).')
