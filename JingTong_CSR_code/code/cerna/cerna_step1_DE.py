# -*- coding: utf-8 -*-
"""
ceRNA step 1: differential expression of mRNA / lncRNA / circRNA in GSE153761
(cervical cartilage endplate, 3 degenerate vs 3 normal) using its native ceRNA array,
then focus DE mRNAs on the Jingtong-Granules hub genes / 89 targets.
"""
import os, io, gzip, numpy as np, pandas as pd
from scipy import stats
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT=os.path.dirname(os.path.abspath(__file__))
ML=os.path.join(ROOT,'ML_GSE223227')

# ---- 1. GPL22120 probe annotation ----
ann={}
hdr=None
with io.open(os.path.join(ROOT,'GPL22120-25936.txt'),encoding='utf-8',errors='replace') as f:
    for l in f:
        if l.startswith(('!','^','#')): continue
        p=l.rstrip('\n').split('\t')
        if hdr is None and p and p[0]=='ID':
            hdr=p; ix={c:i for i,c in enumerate(hdr)}; continue
        if hdr is None: continue
        def g(c): return p[ix[c]] if c in ix and ix[c]<len(p) else ''
        ann[p[0]]=dict(TYPE=g('TYPE'),SYM=g('GENE_SYMBOL'),ACC=g('ACCESSION'),
                       CIRC=g('CIRCRNA_ID'),ASSOC=g('ASSOCIATED_GENE'),ENTREZ=g('ENTREZGENEID'))
print('annotated probes:',len(ann))

# mRNA GENE_SYMBOL is empty -> map ENTREZ -> symbol via mygene
ez=sorted({ann[pr]['ENTREZ'].split('.')[0] for pr in ann if ann[pr]['TYPE']=='mRNA' and ann[pr]['ENTREZ']})
import mygene
mg=mygene.MyGeneInfo()
e2s={}
for r in mg.querymany(ez,scopes='entrezgene',fields='symbol',species='human',verbose=False):
    if 'symbol' in r: e2s[str(r['query'])]=r['symbol']
for pr in ann:
    if ann[pr]['TYPE']=='mRNA' and not ann[pr]['SYM']:
        ann[pr]['SYM']=e2s.get(ann[pr]['ENTREZ'].split('.')[0],'')
print('mapped mRNA symbols:',len(e2s))

# ---- 2. GSE153761 expression matrix (series matrix table) ----
txt=gzip.open(os.path.join(ML,'GSE153761_series_matrix.txt.gz'),'rt',encoding='utf-8',errors='replace').read().splitlines()
beg=[i for i,l in enumerate(txt) if l.startswith('!series_matrix_table_begin')][0]
end=[i for i,l in enumerate(txt) if l.startswith('!series_matrix_table_end')][0]
rows=[l.split('\t') for l in txt[beg+1:end]]
cols=[c.strip('"') for c in rows[0]]
dat=pd.DataFrame(rows[1:],columns=cols).set_index(cols[0])
dat=dat.apply(pd.to_numeric,errors='coerce').dropna()
dat.index=[i.strip('"') for i in dat.index]
samp=list(dat.columns)            # GSM4653870..875 ; first 3 degenerate, last 3 normal
deg,nor=samp[:3],samp[3:]
print('expression matrix:',dat.shape,'| degenerate',deg,'normal',nor)

# ---- 3. DE per RNA type ----
def de(df):
    a=df[deg].values; b=df[nor].values
    lfc=a.mean(1)-b.mean(1)
    p=np.array([stats.ttest_ind(a[i],b[i]).pvalue for i in range(len(df))])
    return pd.DataFrame({'log2FC':lfc,'p':p},index=df.index)

res={}
for typ in ['mRNA','lncRNA','circRNA']:
    probes=[pr for pr in dat.index if ann.get(pr,{}).get('TYPE')==typ]
    sub=dat.loc[probes]
    d=de(sub)
    d['TYPE']=typ
    d['SYM']=[ann[pr]['SYM'] or ann[pr]['ASSOC'] for pr in d.index]
    d['ACC']=[ann[pr]['ACC'] for pr in d.index]
    d['CIRC']=[ann[pr]['CIRC'] for pr in d.index]
    sig=d[(d.p<0.05)&(d.log2FC.abs()>=1)]
    res[typ]=d
    print(f'{typ}: {len(sub)} probes, DE(|log2FC|>=1 & p<0.05)= {len(sig)} (up {sum(sig.log2FC>0)}, down {sum(sig.log2FC<0)})')
    sig.sort_values('p').to_csv(os.path.join(OUT,f'ceRNA_DE_{typ}.csv'),encoding='utf-8-sig')

# ---- 4. focus DE mRNAs on hub genes / 89 targets ----
hub=['AKT1','IL1B','IL6','JUN','EGFR','TP53','MMP9','IL4','IL1A','MMP3','VCAM1']
mR=res['mRNA']
mR_hub=mR[mR.SYM.isin(hub)]
print('\nHub-gene mRNA probes on array (DE direction):')
print(mR_hub[['SYM','log2FC','p']].sort_values('p').to_string())
mR_hub.to_csv(os.path.join(OUT,'ceRNA_hub_mRNA.csv'),encoding='utf-8-sig')
# DE mRNAs (relaxed p<0.05) that are hub genes -> ceRNA mRNA anchors
anchors=sorted(set(mR_hub[mR_hub.p<0.05].SYM)|set(mR_hub[mR_hub.log2FC.abs()>=0.5].SYM))
print('\nceRNA mRNA anchors (hub genes with signal):',anchors)
open(os.path.join(OUT,'ceRNA_anchors.txt'),'w').write(','.join(anchors))
print('\nSaved ceRNA_DE_*.csv, ceRNA_hub_mRNA.csv, ceRNA_anchors.txt')
