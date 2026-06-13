# -*- coding: utf-8 -*-
"""
ADMET / drug-likeness of Jingtong Granules key active compounds.
- fetches canonical SMILES from PubChem by name (accurate, not hand-typed)
- RDKit physicochemical descriptors + Lipinski / Veber / Egan drug-likeness rules
- writes a SwissADME/ADMETlab-ready SMILES file for full PK/tox endpoints
"""
import os, time, urllib.request, urllib.parse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns

OUT=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(OUT,'figures')
COMPOUNDS=['quercetin','kaempferol','stigmasterol','catechin','ginsenoside Rh2',
           'tetrahydropalmatine','tetrahydroberberine','cavidine','cryptopine',
           'ligustilide','puerarin','paeoniflorin']

def pubchem_smiles(name):
    base='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/%s/TXT'
    for prop in ('IsomericSMILES','CanonicalSMILES','SMILES'):
        try:
            url=base%(urllib.parse.quote(name),prop)
            s=urllib.request.urlopen(url,timeout=30).read().decode().strip().splitlines()[0].strip()
            if s and s!='0': return s
        except Exception: time.sleep(0.5); continue
    return None

rows=[]; smi_lines=[]
for name in COMPOUNDS:
    smi=pubchem_smiles(name)
    if not smi: print('no SMILES:',name); continue
    m=Chem.MolFromSmiles(smi)
    if m is None: print('bad SMILES:',name); continue
    mw=Descriptors.MolWt(m); logp=Crippen.MolLogP(m)
    hbd=Lipinski.NumHDonors(m); hba=Lipinski.NumHAcceptors(m)
    tpsa=rdMolDescriptors.CalcTPSA(m); rot=Descriptors.NumRotatableBonds(m)
    arom=rdMolDescriptors.CalcNumAromaticRings(m); fsp3=rdMolDescriptors.CalcFractionCSP3(m)
    lip_v=sum([mw>500,logp>5,hbd>5,hba>10])                 # Lipinski violations
    veber=(rot<=10 and tpsa<=140)
    egan=(tpsa<=131.6 and logp<=5.88)
    gi='High' if (tpsa<=131.6 and -0.4<=logp<=5.88) else 'Low'  # BOILED-Egg-style approx
    rows.append({'compound':name,'MW':round(mw,1),'cLogP':round(logp,2),'HBD':hbd,'HBA':hba,
                 'TPSA':round(tpsa,1),'RotB':rot,'AromRings':arom,'FractionCsp3':round(fsp3,2),
                 'Lipinski_violations':lip_v,'Veber_pass':veber,'Egan_pass':egan,'GI_absorption~':gi})
    smi_lines.append(f'{smi} {name}')
    time.sleep(0.3)

df=pd.DataFrame(rows)
df.to_csv(os.path.join(OUT,'ADMET_druglikeness.csv'),index=False,encoding='utf-8-sig')
open(os.path.join(OUT,'compounds_for_SwissADME.smi'),'w').write('\n'.join(smi_lines))
print(df.to_string(index=False))
print('\nSaved ADMET_druglikeness.csv + compounds_for_SwissADME.smi (submit to SwissADME/ADMETlab for full PK/tox)')

# radar-ish summary figure: Lipinski compliance
# NOTE: the standalone Fig_ADMET heatmap has been RETIRED — its physicochemical descriptors
# (MW/cLogP/TPSA/Rot.bonds/HBD/HBA) are now consolidated into panel (a) of Fig_ADMET_SwissADME
# (see swissadme_figures.py). This script now only produces ADMET_druglikeness.csv +
# compounds_for_SwissADME.smi; no figure is generated here.
print('ADMET table + SMILES written; figure consolidated into Fig_ADMET_SwissADME panel (a).')
