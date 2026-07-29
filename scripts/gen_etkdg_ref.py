#!/usr/bin/env python3
"""RDKit ETKDG-embedded-conformer MMFF energy reference (multi-seed).

Reads scripts/val_set/*.sdf, re-embeds each with ETKDGv3 at K seeds, and records
the MMFF94s single-point energy of each embedded conformer plus the mean/stddev
across seeds. The mean is the strain/quality baseline used by validate_etkdg.py.

Using multiple seeds (default 42,43,100,7 via env WEBMM_SEEDS) removes the
single-seed noise floor: RDKit-vs-RDKit single-seed r~0.99, multi-seed-mean
r~0.997 (see rdkit_self_consistency.py). r=1.0 is NOT the target on this harness.

Writes scripts/etkdg_val/rdkit_ref.json:
  {name: {mean, stddev, seeds:{"42":e,...}, n_embedded, n_atoms, embed_ok}}
"""
import os, glob, json, math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import rdkit
print(f"# rdkit {rdkit.__version__}", flush=True)

IN = 'scripts/val_set'
OUT = 'scripts/etkdg_val'
SEEDS = [int(s) for s in os.environ.get('WEBMM_SEEDS', '42,43,100,7').split(',')]
print(f"# seeds {SEEDS}", flush=True)
os.makedirs(OUT, exist_ok=True)


def embed_energy(mol, seed):
    """Embed (fresh) with ETKDGv3 at `seed`; return MMFF94s energy or None."""
    mol.RemoveAllConformers()
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.useBasicKnowledge = True
    p.useMacrocycleTorsions = True
    try:
        cid = AllChem.EmbedMolecule(mol, p)
    except Exception:
        cid = -1
    if cid < 0:
        return None
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    if mp is None:
        return None
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=0)
    return round(ff.CalcEnergy(), 4)


ref = {}
for f in sorted(glob.glob(os.path.join(IN, '*.sdf'))):
    name = os.path.basename(f).replace('.sdf', '')
    suppl = Chem.SDMolSupplier(f, removeHs=False, sanitize=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        ref[name] = {'error': 'parse'}
        continue
    mol = mols[0]
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
    seeds = {}
    for s in SEEDS:
        e = embed_energy(mol, s)
        if e is not None:
            seeds[str(s)] = e
    n_emb = len(seeds)
    if n_emb == 0:
        ref[name] = {'mean': None, 'stddev': None, 'seeds': {},
                     'n_embedded': 0, 'n_atoms': n, 'embed_ok': False}
        continue
    vals = list(seeds.values())
    mean = sum(vals) / n_emb
    stddev = math.sqrt(sum((v - mean) ** 2 for v in vals) / n_emb)
    ref[name] = {
        'mean': round(mean, 4),
        'stddev': round(stddev, 4),
        'seeds': seeds,
        'n_embedded': n_emb,
        'n_atoms': n,
        'embed_ok': n_emb == len(SEEDS),
    }

json.dump(ref, open(os.path.join(OUT, 'rdkit_ref.json'), 'w'), indent=1)
ok = sum(1 for v in ref.values() if isinstance(v, dict) and v.get('embed_ok'))
print(f"# wrote {OUT}/rdkit_ref.json: {ok}/{len(ref)} fully embedded "
      f"({len(SEEDS)} seeds each)", flush=True)
