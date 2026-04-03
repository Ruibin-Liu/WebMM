#!/usr/bin/env python3
"""
Comprehensive comparison between WebMM and RDKit MMFF94s implementation.
Tests atom typing, energy calculation, and optimization on identical inputs.
"""

from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers, rdDistGeom
import math
import json
import subprocess
import sys

# Test molecules with SDF content (same as used in WebMM tests)
TEST_MOLECULES = {
    "water": """Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END""",
    "methane": """Methane
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6320    0.6320    0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6320   -0.6320   -0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6320    0.6320   -0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6320   -0.6320    0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END""",
    "formaldehyde": """Formaldehyde
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END""",
    "ethanol": """Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END""",
    "acetic_acid": """Acetic acid
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
    -0.9335   -0.0601   -0.2304 C   0  0  0  0  0  0  0  0  0  0  0  0
     0.4936    0.2789    0.0469 C   0  0  0  0  0  0  0  0  0  0  0  0
     1.0325    1.3566   -0.1361 O   0  0  0  0  0  0  0  0  0  0  0  0
     1.1814   -0.7645    0.5462 O   0  0  0  0  0  0  0  0  0  0  0  0
    -1.4427    0.8203   -0.6327 H   0  0  0  0  0  0  0  0  0  0  0  0
    -1.4305   -0.3576    0.6963 H   0  0  0  0  0  0  0  0  0  0  0  0
    -0.9867   -0.8625   -0.9702 H   0  0  0  0  0  0  0  0  0  0  0  0
     2.0859   -0.4112    0.6800 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END""",
    "ethane": """Ethane
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
    -0.7558    0.0071   -0.0165 C   0  0  0  0  0  0  0  0  0  0  0  0
     0.7558   -0.0071    0.0165 C   0  0  0  0  0  0  0  0  0  0  0  0
    -1.1634   -0.1004    0.9931 H   0  0  0  0  0  0  0  0  0  0  0  0
    -1.1223    0.9481   -0.4375 H   0  0  0  0  0  0  0  0  0  0  0  0
    -1.1346   -0.8156   -0.6303 H   0  0  0  0  0  0  0  0  0  0  0  0
     1.1346    0.8156    0.6303 H   0  0  0  0  0  0  0  0  0  0  0  0
     1.1634    0.1004   -0.9931 H   0  0  0  0  0  0  0  0  0  0  0  0
     1.1223   -0.9481    0.4375 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  2  8  1  0
M  END""",
    "ammonia": """Ammonia
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.9500    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4750    0.8227    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4750   -0.8227    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
M  END""",
}


def calc_angle(p1, p2, p3):
    """Calculate angle between three points in degrees."""
    v1 = (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])
    v2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
    dot = sum(a * b for a, b in zip(v1, v2))
    n1 = sum(a**2 for a in v1) ** 0.5
    n2 = sum(a**2 for a in v2) ** 0.5
    if n1 > 0 and n2 > 0:
        cos_theta = max(-1, min(1, dot / (n1 * n2)))
        return math.degrees(math.acos(cos_theta))
    return 0.0


def dist3(p1, p2):
    """Calculate distance between two points."""
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2) ** 0.5


def run_rdkit_test(name, sdf_content):
    """Run RDKit MMFF94s on the same SDF content."""
    mol = Chem.MolFromMolBlock(sdf_content, removeHs=False)
    if mol is None:
        return {"error": "Failed to parse SDF"}

    # Get MMFF94s properties
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")

    # Atom types
    atom_types = []
    for i in range(mol.GetNumAtoms()):
        at = mp.GetMMFFAtomType(i)
        atom_types.append(at)

    # Energy on initial coordinates
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mp, confId=0)
    initial_energy = ff.CalcEnergy()

    # Optimize
    mol_copy = Chem.Mol(mol)
    Chem.AssignStereochemistry(
        mol_copy, cleanIt=True, force=True, flagPossibleStereoCenters=True
    )
    rdForceFieldHelpers.MMFFOptimizeMolecule(
        mol_copy, mmffVariant="MMFF94s", maxIters=500
    )

    # Final energy
    mp2 = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_copy, mmffVariant="MMFF94s")
    ff2 = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_copy, mp2, confId=0)
    final_energy = ff2.CalcEnergy()

    # Get optimized coordinates
    conf = mol_copy.GetConformer()
    coords = []
    for i in range(mol_copy.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])

    # Bond lengths
    bond_lengths = []
    for bond in mol_copy.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        p1 = conf.GetAtomPosition(i)
        p2 = conf.GetAtomPosition(j)
        d = ((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2) ** 0.5
        bond_lengths.append({"i": i, "j": j, "length": round(d, 4)})

    # Angles
    angles = []
    for i in range(mol_copy.GetNumAtoms()):
        for j in range(mol_copy.GetNumAtoms()):
            for k in range(mol_copy.GetNumAtoms()):
                if i != j and j != k and i != k:
                    if mol_copy.GetBondBetweenAtoms(
                        i, j
                    ) and mol_copy.GetBondBetweenAtoms(j, k):
                        p1 = conf.GetAtomPosition(i)
                        p2 = conf.GetAtomPosition(j)
                        p3 = conf.GetAtomPosition(k)
                        angle = calc_angle(
                            [p1.x, p1.y, p1.z], [p2.x, p2.y, p2.z], [p3.x, p3.y, p3.z]
                        )
                        angles.append(
                            {"i": i, "j": j, "k": k, "angle": round(angle, 2)}
                        )

    return {
        "n_atoms": mol.GetNumAtoms(),
        "atom_types": atom_types,
        "initial_energy": round(initial_energy, 6),
        "final_energy": round(final_energy, 6),
        "optimized_coords": coords,
        "bond_lengths": bond_lengths,
        "angles": angles,
    }


def main():
    print("=" * 80)
    print("RDKit MMFF94s Reference Data (from same SDF inputs as WebMM)")
    print("=" * 80)

    all_results = {}

    for name, sdf in TEST_MOLECULES.items():
        print(f"\n{'=' * 40}")
        print(f"  {name.upper()}")
        print(f"{'=' * 40}")

        result = run_rdkit_test(name, sdf)
        all_results[name] = result

        if "error" in result:
            print(f"  ERROR: {result['error']}")
            continue

        print(f"  Atoms: {result['n_atoms']}")
        print(f"  Atom types: {result['atom_types']}")
        print(f"  Initial energy: {result['initial_energy']:.6f} kcal/mol")
        print(f"  Final energy: {result['final_energy']:.6f} kcal/mol")

        print(f"\n  Optimized coordinates:")
        for i, c in enumerate(result["optimized_coords"]):
            print(f"    {i}: {c[0]:.4f} {c[1]:.4f} {c[2]:.4f}")

        print(f"\n  Bond lengths:")
        for b in result["bond_lengths"]:
            print(f"    {b['i']}-{b['j']}: {b['length']:.4f}")

        print(f"\n  Angles:")
        for a in result["angles"]:
            print(f"    {a['i']}-{a['j']}-{a['k']}: {a['angle']:.2f}")

    # Save as JSON for comparison
    with open("/tmp/rdkit_sdf_reference.json", "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"\n\n{'=' * 80}")
    print("Reference data saved to /tmp/rdkit_sdf_reference.json")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
