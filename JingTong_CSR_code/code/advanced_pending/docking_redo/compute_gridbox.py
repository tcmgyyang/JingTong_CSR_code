# -*- coding: utf-8 -*-
"""
Compute AutoDock Vina grid boxes from the verified receptor structures, centred on the
co-crystallised ligand (true binding pocket) where available. Writes vina config_*.txt
and extracts the reference ligand for a re-docking (RMSD<2 A) positive control.
"""
import os, statistics as st
ROOT=os.path.dirname(os.path.abspath(__file__))
AA={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET',
    'PHE','PRO','SER','THR','TRP','TYR','VAL'}
WATER={'HOH','WAT'}
# target -> (pdb, cocrystal ligand resname or None, note)
TARGETS={
 'AKT1':('4EKL','0RF','ATP-competitive pocket (GDC-0068)'),
 'IL1B':('8RYS',None,'IL-1beta (only sulfate bound) -> blind dock over cytokine'),
 'MMP9':('1GKC','NFH','catalytic-Zn / S1prime pocket (reverse-hydroxamate)'),
 'IL6' :('1ALU',None,'no pocket ligand -> blind dock over 4-helix bundle (Site II/III)'),
 'TP53':('2OCJ',None,'no pocket ligand -> blind dock over core/DNA-binding domain (hard target)'),
}
def atoms(path, chain='A'):
    P=[]
    for ln in open(path,encoding='utf-8',errors='replace'):
        if ln[:6].strip() not in ('ATOM','HETATM'): continue
        ch=ln[21]
        if ch not in (chain,' '): continue
        try: x,y,z=float(ln[30:38]),float(ln[38:46]),float(ln[46:54])
        except ValueError: continue
        P.append((ln[:6].strip(),ln[17:20].strip(),x,y,z,ln))
    return P
def centroid(pts):
    return (st.mean(p[2] for p in pts), st.mean(p[3] for p in pts), st.mean(p[4] for p in pts))
def span(pts):
    return (max(p[2] for p in pts)-min(p[2] for p in pts),
            max(p[3] for p in pts)-min(p[3] for p in pts),
            max(p[4] for p in pts)-min(p[4] for p in pts))

print(f'{"target":6s} {"pdb":5s} {"center (x,y,z)":28s} {"box":14s} mode')
for tgt,(pdb,lig,note) in TARGETS.items():
    P=atoms(os.path.join(ROOT,pdb+'.pdb'))
    if lig:
        L=[p for p in P if p[0]=='HETATM' and p[1]==lig]
        if not L:  # try any chain
            L=[p for ln in [pdb] for p in atoms(os.path.join(ROOT,pdb+'.pdb'),chain=' ')+atoms(os.path.join(ROOT,pdb+'.pdb'),'B') if p[0]=='HETATM' and p[1]==lig]
        cx,cy,cz=centroid(L); sx=sy=sz=22.5; mode=f'pocket({lig},{len(L)} atoms)'
        # save reference ligand for redocking RMSD control
        with open(os.path.join(ROOT,f'{tgt}_ref_ligand.pdb'),'w') as f:
            for p in L: f.write(p[5])
    else:
        prot=[p for p in P if p[1] in AA]
        cx,cy,cz=centroid(prot); spx,spy,spz=span(prot)
        sx,sy,sz=min(spx+8,30),min(spy+8,30),min(spz+8,30); mode='blind(whole domain)'
    cfg=(f"# Vina config — {tgt} ({pdb}) — {note}\n"
         f"receptor = {tgt}_receptor.pdbqt\n"
         f"ligand   = quercetin.pdbqt   # also run with kaempferol.pdbqt\n"
         f"center_x = {cx:.2f}\ncenter_y = {cy:.2f}\ncenter_z = {cz:.2f}\n"
         f"size_x = {sx:.0f}\nsize_y = {sy:.0f}\nsize_z = {sz:.0f}\n"
         f"exhaustiveness = 32\nnum_modes = 9\nenergy_range = 4\n")
    open(os.path.join(ROOT,f'config_{tgt}.txt'),'w').write(cfg)
    print(f'{tgt:6s} {pdb:5s} ({cx:7.2f},{cy:7.2f},{cz:7.2f})   {sx:.0f}x{sy:.0f}x{sz:.0f}      {mode}')
print('\nWrote config_AKT1/MMP9/IL6/TP53.txt and *_ref_ligand.pdb (for redocking control)')
