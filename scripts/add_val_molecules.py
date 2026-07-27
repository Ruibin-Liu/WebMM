#!/usr/bin/env python3
"""Add new validation molecules targeting gaps the energy-debugging exposed.
Generates only NEW molecules, appends SDFs + ref entries to the existing set."""
import json, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# new SMILES targeting: sp2-3-sub (OOP), X-OH (DFSB), fused rings (sb_type),
# halogens, Si, strained rings, zwitterions, larger drugs
NEW = {
    # sp2 centers with 3 substituents (stress OOP cyclic terms)
    'formamidine': 'NC=N', 'acetamidine': 'CC(=N)N', 'biguanide': 'NC(=N)NC(=N)N',
    # X-OH (stress DFSB asymmetric-row stretch-bend)
    'phosphoric_acid': 'O=P(O)(O)O', 'sulfuric_acid': 'OS(=O)(=O)O',
    'p_toluene_sulfonic_acid': 'Cc1ccc(S(=O)(=O)O)cc1',
    # fused heterocycles (stress sb_type canonicalization)
    'benzimidazole': 'c1ccc2[nH]cnc2c1', 'benzothiazole': 'c1ccc2scnc2c1',
    'xanthine': 'O=c1[nH]c(=O)c2[nH]cnc2[nH]1',
    'hypoxanthine': 'O=c1[nH]cnc2c[nH]c12',
    # halogen diversity (stress C-X bond params)
    'trifluoromethane': 'FC(F)F', 'dibromomethane': 'BrCBr', 'diiodomethane': 'ICI',
    'vinyl_chloride': 'C=CCl', 'acetyl_chloride': 'CC(=O)Cl',
    # Si (new element)
    'tetramethylsilane': 'C[Si](C)(C)C',
    # strained / exotic geometry
    'cyclopropene': 'C1=CC1', 'allene': 'C=C=C',
    'ethylene_oxide': 'C1CO1', 'aziridine': 'C1CN1',
    # zwitterion / charge
    'glycine_zwitterion': '[NH3+]CC(=O)[O-]',
    # larger drug-like
    'ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
}

REF = 'scripts/val_set/rdkit_ref.json'
ref = json.load(open(REF))
added, skipped = [], []

for name, smi in NEW.items():
    if name in ref:
        skipped.append(name + '(exists)')
        continue
    try:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        if m is None:
            skipped.append(name + '(parse)'); continue
        p = AllChem.ETKDGv3(); p.randomSeed = 42; p.maxIterations = 200
        if AllChem.EmbedMolecule(m, p) != 0:
            skipped.append(name + '(embed)'); continue
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        Chem.SetAromaticity(m)
        mp = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
        if mp is None:
            skipped.append(name + '(no-MMFF)'); continue
        ff = AllChem.MMFFGetMoleculeForceField(m, mp, confId=0)
        open('scripts/val_set/%s.sdf' % name, 'w').write(Chem.MolToMolBlock(m) + '$$$$\n')
        ref[name] = {
            'source': 'curated',
            'types': [mp.GetMMFFAtomType(i) for i in range(m.GetNumAtoms())],
            'charges': [round(mp.GetMMFFPartialCharge(i), 4) for i in range(m.GetNumAtoms())],
            'energy': round(ff.CalcEnergy(), 4),
            'n_atoms': m.GetNumAtoms(),
        }
        added.append(name)
    except Exception as e:
        skipped.append('%s(%s)' % (name, str(e)[:30]))

json.dump(ref, open(REF, 'w'))
print('ADDED (%d): %s' % (len(added), ', '.join(added)))
print('SKIPPED (%d): %s' % (len(skipped), ', '.join(skipped)))
print('total molecules now:', len(ref))
