"""Extract MMFFANG, MMFFSB, MMFFDef tables from RDKit."""

import sys

sys.path.insert(0, "/opt/homebrew/lib/python3.14/site-packages")

from rdkit import Chem
from rdkit.Chem import AllChem

# Access the internal parameter tables
from rdkit.ForceField import MMFF


# Get default parameter collections
def extract_tables():
    # MMFFProp table - has eqLevel, crd, val, mltb, etc.
    # MMFFDef table - has eqLevel[4] per atom type
    # MMFFAngle table - angle parameters
    # MMFFStbn table - stretch-bend parameters

    # Let's use the Python-level access to the parameter tables
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    mp = AllChem.MMFFGetMoleculeProperties(mol, "MMFF94s")

    # Print MMFFProp entries for atom types 1-99
    print("MMFFProp table (atom types with crd, val, mltb, arom, sbmb):")
    for i in range(1, 100):
        try:
            prop = mp.GetMMFFAtomType(i)
            if prop is not None:
                print(f"  type {i}: {prop}")
        except:
            pass

    # Access the raw parameter tables
    print("\nTrying to access parameter collections...")

    # Look at what's available
    mmff_module = dir(MMFF)
    print(
        f"MMFF module attributes: {[x for x in mmff_module if not x.startswith('_')]}"
    )


extract_tables()
