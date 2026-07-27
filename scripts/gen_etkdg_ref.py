#!/usr/bin/env python3
"""RDKit ETKDG-embedded-conformer MMFF energy reference (mirrors MMFF val-set shape).

Reads scripts/val_set/*.sdf (connectivity source — WebMM re-embeds anyway), re-embeds
each with ETKDGv3 (randomSeed=42, useMacrocycleTorsions), and records the MMFF94s
single-point energy of the embedded conformer — the strain/quality baseline.

Writes scripts/etkdg_val/rdkit_ref.json {name: {energy, n_atoms, embed_ok}}.
"""
import os, glob, json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import rdkit
print(f"# rdkit {rdkit.__version__}", flush=True)

IN = 'scripts/val_set'
OUT = 'scripts/etkdg_val'
os.makedirs(OUT, exist_ok=True)

ref = {}
for f in sorted(glob.glob(os.path.join(IN, '*.sdf'))):
    name = os.path.basename(f).replace('.sdf', '')
    suppl = Chem.SDMolSupplier(f, removeHs=False, sanitize=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        ref[name] = {'error': 'parse'}
        continue
    mol = mols[0]
    # Full sanitize (with kekulize) so aromaticity is perceived from the kekulé
    # SDF bonds; fall back to skip-kekulize for aromatics that can't kekulize.
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            Chem.SetAromaticity(mol)
        except Exception as e:
            ref[name] = {'error': f'sanitize:{e}'}
            continue
    try:
        Chem.FastFindRings(mol)
    except Exception:
        pass
    n = mol.GetNumAtoms()
    # Fresh ETKDGv3 embed (drop any existing conformer from the SDF).
    try:
        mol.RemoveAllConformers()
        p = AllChem.ETKDGv3()
        p.randomSeed = 42
        p.useBasicKnowledge = True
        p.useMacrocycleTorsions = True
        cid = AllChem.EmbedMolecule(mol, p)
        embed_ok = cid >= 0
    except Exception:
        embed_ok = False
    if not embed_ok:
        ref[name] = {'energy': None, 'n_atoms': n, 'embed_ok': False}
        continue
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    if mp is None:
        ref[name] = {'energy': None, 'n_atoms': n, 'embed_ok': False}
        continue
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=0)
    ref[name] = {'energy': round(ff.CalcEnergy(), 4), 'n_atoms': n, 'embed_ok': True}

json.dump(ref, open(os.path.join(OUT, 'rdkit_ref.json'), 'w'), indent=1)
ok = sum(1 for v in ref.values() if isinstance(v, dict) and v.get('embed_ok'))
print(f"# wrote {OUT}/rdkit_ref.json: {ok}/{len(ref)} embedded", flush=True)
