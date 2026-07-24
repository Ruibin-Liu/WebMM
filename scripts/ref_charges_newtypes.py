"""Extract RDKit MMFF94s reference atom types and partial charges for WebMM tests."""
from rdkit import Chem
from rdkit.Chem import AllChem

MOLS = [
    ("dmso", "CS(=O)C"),
    ("sulfone", "CS(=O)(=O)C"),
    ("nitromethane", "C[N+](=O)[O-]"),
    ("benzaldehyde", "O=Cc1ccccc1"),
    ("sulfonamide", "CS(=O)(=O)N"),
    ("acetone", "CC(=O)C"),
]

for name, smi in MOLS:
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    p = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    print(name, [p.GetMMFFAtomType(i) for i in range(mol.GetNumAtoms())])
    print("  ", [round(p.GetMMFFPartialCharge(i), 4) for i in range(mol.GetNumAtoms())])
