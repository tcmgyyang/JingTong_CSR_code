# -*- coding: utf-8 -*-
"""
STITCH supplement: high-confidence chemical-protein targets for the 7 Jingtong-Granules
herbs' active compounds. Queries STITCH (combined score >= 700), maps Ensembl-protein
partners to gene symbols, aggregates to herb level. One of FOUR target sources
(TCMSP + SwissTargetPrediction + HERB + STITCH) that will be merged.
"""
import os, time, urllib.request, urllib.parse
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.dirname(os.path.abspath(__file__))

d = pd.read_csv(os.path.join(ROOT, 'allTargets.symbol.xls'), sep='\t')
seven = ['baishao', 'chuanxiong', 'gegen', 'qianghuo', 'sanqi', 'weilingxian', 'yanhusuo']
sub = d[d.Drug.isin(seven)][['Drug', 'MolId', 'MolName']].drop_duplicates('MolId')
herbmap = sub.groupby('MolName')['Drug'].apply(lambda s: sorted(set(s))).to_dict()
compounds = sorted(set(sub['MolName'].astype(str)))
print('querying STITCH for', len(compounds), 'compounds (score>=700)...')

SCORE = 700
UA = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) jingtong-netpharm/1.0'
def stitch_partners(name, retries=3):
    u = ('http://stitch.embl.de/api/tsv/interactors?identifier=%s&species=9606'
         '&required_score=%d&limit=500&caller_identity=jingtong_netpharm'
         % (urllib.parse.quote(name), SCORE))
    for k in range(retries):
        try:
            req = urllib.request.Request(u, headers={'User-Agent': UA})
            txt = urllib.request.urlopen(req, timeout=45).read().decode('utf-8', 'replace')
            ensp = [ln.strip().split('.')[-1] for ln in txt.splitlines() if 'ENSP' in ln]
            return sorted(set(ensp))
        except Exception:
            time.sleep(5.0 * (k + 1))      # back off progressively
    return None

rows = []          # compound, herb, ENSP
miss = []
for c in compounds:
    p = stitch_partners(c)
    if p is None:
        miss.append(c); print('  502/err:', c); time.sleep(1.0); continue
    for ensp in p:
        for h in herbmap.get(c, []):
            rows.append([c, h, ensp])
    print(f'  {c[:40]:40s} {len(p)} partners')
    time.sleep(3.5)                          # be polite to STITCH (avoid rate-limit/block)

ce = pd.DataFrame(rows, columns=['MolName', 'herb', 'ENSP']).drop_duplicates()
print('\nraw ENSP edges:', len(ce), '| unique ENSP:', ce['ENSP'].nunique(), '| missed compounds:', len(miss))

# map Ensembl protein -> gene symbol
import mygene
mg = mygene.MyGeneInfo()
ensps = sorted(ce['ENSP'].unique())
res = mg.querymany(ensps, scopes='ensembl.protein', fields='symbol', species='human', verbose=False)
e2s = {r['query']: r['symbol'] for r in res if 'symbol' in r}
ce['gene'] = ce['ENSP'].map(e2s)
ce = ce.dropna(subset=['gene'])
if len(ce) == 0:
    open(os.path.join(OUT, 'STITCH_targets.csv'), 'w', encoding='utf-8-sig').write('herb,gene,MolName\n')
    raise SystemExit('STITCH returned no partners (server unavailable) -> empty STITCH_targets.csv')
ce['gene'] = ce['gene'].astype(str).str.upper()
ce[['herb', 'gene', 'MolName']].drop_duplicates().to_csv(os.path.join(OUT, 'STITCH_targets.csv'), index=False, encoding='utf-8-sig')
print('STITCH unique gene targets (7 herbs):', ce['gene'].nunique())
print('per herb:'); print(ce.groupby('herb')['gene'].nunique().to_string())
open(os.path.join(OUT, 'STITCH_missed.txt'), 'w', encoding='utf-8').write('\n'.join(miss))
print('saved STITCH_targets.csv')
