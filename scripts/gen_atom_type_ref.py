#!/usr/bin/env python3
"""Generate RDKit MMFF94s atom types for all test molecules -> JSON reference."""
import json, glob, os, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

out = {}
dirs = ['test_molecules', 'scripts/test_molecules']
files = []
for d in dirs:
    files.extend(glob.glob(os.path.join(d, '*.sdf')))

for f in sorted(set(files)):
    name = os.path.basename(f).replace('.sdf', '')
    with open(f) as fh:
        txt = fh.read()
    # SDMolSupplier handles multi-mol; take first
    suppl = Chem.SDMolSupplier(f, removeHs=False, sanitize=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        out[name] = {'error': 'parse failed'}
        continue
    entries = []
    for mol in mols:
        # Proper sanitization is essential: RDKit's MMFF typing depends on
        # aromaticity perception, which only happens during sanitization.
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            try:
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                Chem.SetAromaticity(mol)
            except Exception:
                pass
        try:
            Chem.FastFindRings(mol)
        except Exception:
            pass
        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        if mp is None:
            entries.append({'error': 'no MMFF props'})
            continue
        types = [mp.GetMMFFAtomType(i) for i in range(mol.GetNumAtoms())]
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        charges = [a.GetFormalCharge() for a in mol.GetAtoms()]
        bonds = []
        for b in mol.GetBonds():
            bt = b.GetBondType()
            order = {Chem.BondType.SINGLE:1, Chem.BondType.DOUBLE:2,
                     Chem.BondType.TRIPLE:3, Chem.BondType.AROMATIC:4}.get(bt, 1)
            bonds.append([b.GetBeginAtomIdx(), b.GetEndAtomIdx(), order])
        # RDKit aromatic atoms
        arom = [int(a.GetIsAromatic()) for a in mol.GetAtoms()]
        entries.append({'types': types, 'syms': syms, 'charges': charges,
                        'bonds': bonds, 'aromatic': arom})
    out[name] = entries

with open('scripts/rdkit_atom_types.json', 'w') as fh:
    json.dump(out, fh, indent=1)
total = sum(len(v) for v in out.values() if isinstance(v, list))
print(f'wrote scripts/rdkit_atom_types.json: {len(out)} files, {total} mols')
