"""Compare WebMM vs RDKit embedded geometry for one molecule (same seed).

Usage: python3 scripts/diag_compare.py <mol_name> [seed]
Prints bond lengths and angles for RDKit's ETKDGv3 embedding (matching
gen_etkdg_ref.py) plus MMFF94s energy; WebMM geometry comes from the
diag_outliers example output (run manually or via subprocess below).
"""
import subprocess, sys, json
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

name = sys.argv[1]
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 42
path = f'scripts/val_set/{name}.sdf'

# --- RDKit embed (mirror gen_etkdg_ref.py) ---
mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
Chem.SanitizeMol(mol)
mol = Chem.AddHs(mol)
p = AllChem.ETKDGv3()
p.randomSeed = seed
p.useRandomCoords = False
ok = AllChem.EmbedMolecule(mol, p)
conf = mol.GetConformer()
ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s'), confId=0)
print(f'# RDKit name={name} seed={seed} embed_ok={ok} mmff={ff.CalcEnergy():.4f}')
conf = mol.GetConformer()
print('# RDKit coords (idx sym x y z)')
for a in mol.GetAtoms():
    p1 = conf.GetAtomPosition(a.GetIdx())
    print(f'C\t{a.GetIdx()}\t{a.GetSymbol()}\t{p1.x:.4f}\t{p1.y:.4f}\t{p1.z:.4f}')
print('# RDKit bonds (i sym_i - j sym_j len)')
for b in mol.GetBonds():
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    p1 = conf.GetAtomPosition(i); p2 = conf.GetAtomPosition(j)
    d = (p1 - p2).Length()
    print(f'B\t{i}\t{mol.GetAtomWithIdx(i).GetSymbol()}\t{j}\t{mol.GetAtomWithIdx(j).GetSymbol()}\t{d:.4f}')
print('# RDKit angles (i-j-k deg)')
for a in mol.GetAtoms():
    j = a.GetIdx()
    nbrs = [n.GetIdx() for n in a.GetNeighbors()]
    for i in nbrs:
        for k in nbrs:
            if i < k:
                ang = rdMolTransforms.GetAngleDeg(conf, i, j, k)
                print(f'ANG\t{i}\t{j}\t{k}\t{ang:.2f}')

# --- WebMM embed via diag_outliers ---
r = subprocess.run(['./target/release/examples/diag_outliers', name, str(seed)],
                   capture_output=True, text=True)
print('--- WebMM ---')
print(r.stdout, end='')
if r.stderr:
    print(r.stderr[:500])
