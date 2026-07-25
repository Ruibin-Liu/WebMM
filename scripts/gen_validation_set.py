#!/usr/bin/env python3
"""Generate a large validation set: curated SMILES + PubChem drugs -> 3D SDFs + RDKit reference.

Outputs:
  scripts/val_set/<name>.sdf        one explicit-H 3D molecule each
  scripts/val_set/rdkit_ref.json    {name: {types, charges, energy, n_atoms}}
"""
import os, sys, json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

OUT = 'scripts/val_set'
os.makedirs(OUT, exist_ok=True)

# ---- Curated SMILES covering the MMFF chemistry space ----
SMILES = {
    # alkanes / alkenes / alkynes
    'ethane': 'CC', 'propane': 'CCC', 'isobutane': 'CC(C)C', 'neopentane': 'CC(C)(C)C',
    'ethylene': 'C=C', 'propene': 'C=CC', 'butadiene': 'C=CC=C', 'isoprene': 'C=C(C)C=C',
    'acetylene': 'C#C', 'propyne': 'CC#C',
    # O: alcohols, ethers, phenols
    'methanol': 'CO', 'ethanol': 'CCO', 'phenol': 'c1ccccc1O', 'catechol': 'c1ccccc1(O)O',
    'dimethyl_ether': 'COC', 'diethyl_ether': 'CCOCC', 'anisole': 'c1ccccc1OC',
    'THF': 'C1CCOC1', 'dioxane': 'C1OCCOC1',
    # O: carbonyls
    'formaldehyde': 'C=O', 'acetaldehyde': 'CC=O', 'benzaldehyde': 'c1ccccc1C=O',
    'acetone': 'CC(=O)C', 'cyclohexanone': 'O=C1CCCCC1',
    # O: acids, esters, anhydride
    'formic_acid': 'OC=O', 'acetic_acid': 'CC(O)=O', 'benzoic_acid': 'c1ccccc1C(O)=O',
    'oxalic_acid': 'OC(=O)C(O)=O', 'methyl_formate': 'COC=O', 'ethyl_acetate': 'CC(=O)OCC',
    'methyl_salicylate': 'O=C(OC)c1ccccc1O', 'acetic_anhydride': 'CC(=O)OC(C)=O',
    # N: amines
    'methylamine': 'CN', 'dimethylamine': 'CNC', 'trimethylamine': 'CN(C)C',
    'tetramethylammonium': 'C[N+](C)(C)C', 'aniline': 'c1ccccc1N',
    'N_methylaniline': 'c1ccccc1NC', 'nitrobenzene': 'c1ccccc1[N+](=O)[O-]',
    'nitromethane': 'C[N+](=O)[O-]',
    # N: amides, imines, nitriles
    'formamide': 'NC=O', 'acetamide': 'CC(N)=O', 'urea': 'NC(N)=O',
    'DMF': 'CN(C)C=O', 'N_methylacetamide': 'CC(=O)NC', 'benzamide': 'c1ccccc1C(N)=O',
    'acetaldimine': 'CC=N', 'benzaldimine': 'c1ccccc1C=N',
    'acetonitrile': 'CC#N', 'benzonitrile': 'c1ccccc1C#N',
    # N heterocycles
    'pyridine': 'c1ccncc1', 'pyrrole': 'c1cc[nH]c1', 'imidazole': 'c1cnc[nH]1',
    'pyrazole': 'c1cc[nH]n1', 'pyrimidine': 'c1cncnc1', 'pyrazine': 'c1cnccn1',
    'indole': 'c1ccc2[nH]ccc2c1', 'quinoline': 'c1ccc2ncccc2c1', 'isoquinoline': 'c1ccc2cnccc2c1',
    'purine': 'c1nc2c(c1)[nH]cnc2', 'pyridazine': 'c1ccnnn1',
    # S chemistry
    'methanethiol': 'CS', 'dimethyl_sulfide': 'CSC', 'dimethyl_disulfide': 'CSSC',
    'thiophene': 'c1ccsc1', 'thiazole': 'c1cncs1', 'benzothiophene': 'c1ccc2sccc2c1',
    'DMSO': 'CS(=O)C', 'dimethyl_sulfone': 'CS(=O)(=O)C', 'methanesulfonamide': 'CS(=O)(=O)N',
    'methanesulfonic_acid': 'CS(=O)(=O)O', 'thioacetone': 'CC(=S)C', 'carbon_disulfide': 'S=C=S',
    # P chemistry
    'trimethylphosphine': 'CP(C)C', 'trimethylphosphine_oxide': 'CP(=O)(C)C',
    'triethyl_phosphate': 'CCOP(=O)(OCC)OCC', 'methylphosphonic_acid': 'CP(=O)(O)O',
    # halogens
    'fluoromethane': 'CF', 'chloroform': 'ClC(Cl)Cl', 'bromobenzene': 'c1ccccc1Br',
    'iodobenzene': 'c1ccccc1I', 'dichloromethane': 'ClCCl', 'carbon_tetrachloride': 'ClC(Cl)(Cl)Cl',
    # rings
    'cyclopropane': 'C1CC1', 'cyclobutane': 'C1CCC1', 'cyclopentane': 'C1CCCC1',
    'cyclohexane': 'C1CCCCC1', 'cycloheptane': 'C1CCCCCC1',
    'cyclopentene': 'C1=CCCC1', 'cyclohexene': 'C1=CCCCC1', 'norbornane': 'C1CC2CCC1C2',
    'adamantane': 'C1C2CC3CC1CC(C2)C3',
    'benzene': 'c1ccccc1', 'naphthalene': 'c1ccc2ccccc2c1', 'anthracene': 'c1ccc2cc3ccccc3cc2c1',
    'furan': 'c1ccoc1', 'tetrahydrofuran': 'C1CCOC1',
    # ions / charged
    'ammonium': '[NH4+]', 'acetate': 'CC(=O)[O-]', 'benzoate': 'c1ccccc1C(=O)[O-]',
    'guanidinium': 'NC(=N)[NH2+]', 'trimethylammonium': 'C[NH+](C)C',
    'imidazolium': 'c1c[nH+]cc1',
    # fused N-hetero / drug-like
    'caffeine': 'Cn1c(=O)c2c(ncn2C)n(C)c1=O', 'theophylline': 'Cn1c(=O)c2[nH]cnc2n(C)c1=O',
    'adenine': 'Nc1ncnc2[nH]cnc12', 'guanine': 'Nc1[nH]c(=O)c2ncnc12', 'uracil': 'O=c1cc[nH]c(=O)n1',
    'aspirin': 'CC(=O)Oc1ccccc1C(O)=O', 'paracetamol': 'CC(=O)Nc1ccc(O)cc1',
    'nicotine': 'c1ccccc1C1CCN(C)CC1', 'histidine_sidechain': 'NC(Cc1c[nH]cn1)=O',
}

# ---- Real drugs fetched from PubChem by name ----
DRUGS = [
    'ibuprofen', 'naproxen', 'ketoprofen', 'diclofenac', 'warfarin',
    'metformin', 'lidocaine', 'propranolol', 'atenolol', 'captopril',
    'ranitidine', 'famotidine', 'omeprazole', 'diazepam', 'fluoxetine',
    'sertraline', 'amitriptyline', 'chlorpromazine', 'penicillin G',
    'cefalexin', 'succinylsulfathiazole', 'chloramphenicol', 'acyclovir',
    'zidovudine', 'methotrexate', 'ritonavir',
]

def fetch_pubchem_sdf(name):
    """Fetch via curl (subprocess/urllib import hangs in some sandboxes)."""
    import subprocess  # lazy; may hang in restricted environments
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF?record_type=3d'
    try:
        r = subprocess.run(['curl', '-s', '--max-time', '20', url],
                           capture_output=True, text=True, timeout=25)
        return r.stdout if ('V2000' in r.stdout or 'V3000' in r.stdout) else None
    except Exception as e:
        print(f'  PubChem fetch failed for {name}: {e}', file=sys.stderr)
        return None

def process(name, mol, source):
    """Add Hs, 3D embed, sanitize; return SDF + RDKit ref or None."""
    try:
        m = Chem.AddHs(mol)
        parm = AllChem.ETKDGv3(); parm.randomSeed = 42; parm.maxIterations = 200
        if AllChem.EmbedMolecule(m, parm) != 0:
            return None
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        Chem.SetAromaticity(m)
        mp = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
        if mp is None:
            return None
        ff = AllChem.MMFFGetMoleculeForceField(m, mp, confId=0)
        sdf = Chem.MolToMolBlock(m)
        ref = {
            'source': source,
            'types': [mp.GetMMFFAtomType(i) for i in range(m.GetNumAtoms())],
            'charges': [round(mp.GetMMFFPartialCharge(i), 4) for i in range(m.GetNumAtoms())],
            'energy': round(ff.CalcEnergy(), 4),
            'n_atoms': m.GetNumAtoms(),
        }
        return sdf, ref
    except Exception as e:
        print(f'  process failed for {name}: {e}', file=sys.stderr)
        return None

ref = {}
fetch_drugs = '--drugs' in sys.argv
# curated SMILES
n_ok = n_fail = 0
for name, smi in SMILES.items():
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f'  bad SMILES: {name}', file=sys.stderr); n_fail += 1; continue
    r = process(name, mol, 'curated')
    if r:
        sdf, d = r; open(f'{OUT}/{name}.sdf', 'w').write(sdf); ref[name] = d; n_ok += 1
        print(f'  curated {name} ok ({d["n_atoms"]} atoms)', flush=True)
    else:
        n_fail += 1
print(f'-- curated done: {n_ok} ok, {n_fail} fail --', flush=True)
json.dump(ref, open(f'{OUT}/rdkit_ref.json', 'w'), indent=1)

# PubChem drugs (opt-in: slow, network)
n_drug = 0
if fetch_drugs:
    for name in DRUGS:
        sdf_txt = fetch_pubchem_sdf(name)
        if not sdf_txt:
            continue
        mol = Chem.MolFromMolBlock(sdf_txt, sanitize=False)
        if mol is None:
            continue
        r = process('drug_' + name, mol, 'pubchem')
        if r:
            sdf, d = r; open(f'{OUT}/drug_{name}.sdf', 'w').write(sdf); ref['drug_' + name] = d; n_drug += 1
            print(f'  drug {name} ok ({d["n_atoms"]} atoms)', flush=True)
    json.dump(ref, open(f'{OUT}/rdkit_ref.json', 'w'), indent=1)

print(f'curated OK: {n_ok}, fail: {n_fail}, drugs: {n_drug}, total: {len(ref)}')
