#!/bin/bash
# =============================================================================
# build_complex.sh — build a GROMACS protein-ligand complex from a docked pose,
# UNIFORMLY for AKT1/IL6/IL1B (AMBER ff99SB-ILDN protein + GAFF2 ligand via acpype +
# TIP3P water). Produces merge_gmx.gro / merge_topol.top / lig.gro / lig.itp that
# feed directly into MD_scripts/01_run_md.sh.
# Requires: AmberTools (antechamber/parmchk2), acpype, GROMACS, Open Babel.
#   conda install -c conda-forge ambertools acpype gromacs openbabel
# Usage: bash build_complex.sh <TARGET> <receptor.pdb> <ligand_pose.pdbqt|mol2|pdb> <outdir>
# Example: bash build_complex.sh AKT1 4EKL.pdb results/AKT1_quercetin_out.pdbqt 分子动力-AKT1-que
# NOTE: MMP9 (catalytic Zn) and TP53 are docking-only; do NOT build MD complexes for them here.
# =============================================================================
set -euo pipefail
TGT=${1:?target}; REC=$(realpath "${2:?receptor pdb}"); POSE=$(realpath "${3:?ligand pose}"); OUT=${4:?outdir}
FF=${FF:-amber99sb-ildn}; WATER=${WATER:-tip3p}
mkdir -p "$OUT"; cd "$OUT"; G=${GMX:-gmx}

# ---- 1. clean receptor (protein only: drop HETATM/water/ions) ----
grep -E '^ATOM|^TER|^END' "$REC" > receptor_clean.pdb
$G pdb2gmx -f receptor_clean.pdb -o protein.gro -p topol.top -i posre.itp -water $WATER -ff $FF -ignh

# ---- 2. ligand: best docked pose -> mol2 (single model) -> acpype (GAFF2, AM1-BCC) ----
case "$POSE" in
  *.pdbqt) obabel "$POSE" -O lig_pose.mol2 -f 1 -l 1 ;;   # first (best) mode
  *)       obabel "$POSE" -O lig_pose.mol2 ;;
esac
acpype -i lig_pose.mol2 -b MOL -n 0 -a gaff2 -c bcc
cp MOL.acpype/MOL_GMX.itp lig_full.itp
cp MOL.acpype/MOL_GMX.gro lig.gro

# ---- 3. split acpype itp into atomtypes (-> lig_atomtypes.itp) and the rest (-> lig.itp) ----
awk 'BEGIN{at=0}
 /^\[ *atomtypes *\]/{at=1; print > "lig_atomtypes.itp"; next}
 at==1 && /^\[/ && !/atomtypes/{at=2}
 at==1{print > "lig_atomtypes.itp"; next}
 {print > "lig.itp"}' lig_full.itp

# ---- 4. assemble complex coordinates (protein + ligand) ----
pn=$(sed -n '2p' protein.gro); ln=$(sed -n '2p' lig.gro)
tot=$(( pn + ln ))
{ echo "Protein-ligand complex"; echo " $tot";
  tail -n +3 protein.gro | head -n $pn;
  tail -n +3 lig.gro     | head -n $ln;
  tail -n 1 protein.gro; } > merge_gmx.gro

# ---- 5. assemble merged topology ----
cp topol.top merge_topol.top
# 5a. insert ligand atomtypes right after the forcefield include
awk '{print}
 /#include ".*forcefield.itp"/ && !done{print "\n; ligand atomtypes\n#include \"lig_atomtypes.itp\""; done=1}' merge_topol.top > t1 && mv t1 merge_topol.top
# 5b. insert ligand topology + POSRES_LIG block just before [ system ]
awk '/^\[ *system *\]/ && !d{print "; ligand topology\n#include \"lig.itp\"\n#ifdef POSRES_LIG\n#include \"posre_lig.itp\"\n#endif\n"; d=1}1' merge_topol.top > t2 && mv t2 merge_topol.top
# 5c. add ligand to [ molecules ]
echo "MOL                 1" >> merge_topol.top

echo ">>> [$OUT] complex built: merge_gmx.gro + merge_topol.top (+ lig.gro/lig.itp)"
echo ">>> next: bash MD_scripts/01_run_md.sh $OUT <gpu_id>"
