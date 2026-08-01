"""Measure MMFF implementation divergence: webmm-MMFF vs RDKit-MMFF on the SAME geometry.

For each molecule in a list: embed with RDKit ETKDGv3 (reference protocol), dump coords,
score with RDKit MMFF94s AND webmm MMFF94s (via diag_outliers external-coords mode),
and report the divergence |E_webmm - E_rdkit| on the identical geometry.

Usage: python3 scripts/diag_mmff_divergence.py
"""
import json, math, subprocess, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

MOLS = sys.argv[1:] or [
    'phosphoric_acid', 'trimethylphosphine_oxide', 'triethyl_phosphate', 'methylphosphonic_acid',
    'dimethyl_disulfide', 'carbon_disulfide', 'thioacetone', 'dimethyl_sulfone',
    'naphthalene', 'xanthine', 'theophylline', 'caffeine',
    'benzene', 'paracetamol', 'ibuprofen', 'DMF', 'anthracene', 'glycine_zwitterion',
    'trimethylphosphine', 'phosphine', 'triphenylphosphine_oxide',
]

rows = []
for name in MOLS:
    path = f'scripts/val_set/{name}.sdf'
    try:
        mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
        Chem.SanitizeMol(mol)
    except Exception:
        continue
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    p.useBasicKnowledge = True
    p.useMacrocycleTorsions = True
    try:
        cid = AllChem.EmbedMolecule(mol, p)
    except Exception:
        cid = -1
    if cid < 0:
        print(f'{name:<26} EMBED-FAIL')
        continue
    conf = mol.GetConformer()
    # write coords for webmm
    with open('/tmp/mmff_div_coords.txt', 'w') as f:
        for a in mol.GetAtoms():
            pos = conf.GetAtomPosition(a.GetIdx())
            f.write(f'C\t{a.GetIdx()}\t{a.GetSymbol()}\t{pos.x:.6f}\t{pos.y:.6f}\t{pos.z:.6f}\n')
    # RDKit MMFF94s
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    if mp is None:
        print(f'{name:<26} NO-MMFF-PARAMS')
        continue
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=0)
    e_rd = ff.CalcEnergy()
    # webmm MMFF94s
    r = subprocess.run(['./target/release/examples/diag_outliers', name, '42', '/tmp/mmff_div_coords.txt'],
                       capture_output=True, text=True)
    line = [l for l in r.stdout.splitlines() if l.startswith('# webmm-MMFF')]
    if not line:
        print(f'{name:<26} WEBMM-ERR: {r.stdout.strip()[:60]}')
        continue
    e_web = float(line[0].split(':')[-1].split('|')[0].strip())
    rows.append((name, e_rd, e_web, e_web - e_rd))
    print(f'{name:<26} rdkit={e_rd:>10.2f} webmm={e_web:>10.2f} d={e_web-e_rd:>+9.2f}')

print(f'\nmax |d| = {max(abs(r[3]) for r in rows):.2f} kcal (mean {sum(abs(r[3]) for r in rows)/len(rows):.2f})')
