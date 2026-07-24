#!/usr/bin/env python3
"""
Direct comparison of RDKit MMFF94s vs WebMM MMFF94s optimization.

Generates a Rust test file that:
1. Parses the same SDF as RDKit used
2. Uses RDKit's optimized coords as starting point for WebMM
3. Compares WebMM's optimized geometry against RDKit's

Then runs the test and analyzes results.
"""

import json
import subprocess
import sys
import os
import math

from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers, rdDistGeom, AllChem

SMILES_LIST = [
    ("water", "O"),
    ("methane", "C"),
    ("ammonia", "N"),
    ("ethane", "CC"),
    ("ethanol", "CCO"),
    ("formaldehyde", "C=O"),
    ("acetic_acid", "CC(=O)O"),
    ("acetone", "CC(=O)C"),
    ("methylamine", "CN"),
    ("dimethyl_ether", "COC"),
    ("propane", "CCC"),
    ("benzene", "c1ccccc1"),
    ("phenol", "c1ccc(O)cc1"),
    ("acetaldehyde", "CC=O"),
    ("glycine", "NCC(=O)O"),
    ("cyclohexane", "C1CCCCC1"),
]


def get_rdkit_optimized(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if mp is None:
        return None

    atom_types = [mp.GetMMFFAtomType(i) for i in range(mol.GetNumAtoms())]

    mol_opt = Chem.RWMol(mol)
    res = rdForceFieldHelpers.MMFFOptimizeMolecule(
        mol_opt, mmffVariant="MMFF94s", maxIters=2000
    )

    mp2 = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_opt, mmffVariant="MMFF94s")
    ff_opt = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_opt, mp2, confId=0)
    final_energy = ff_opt.CalcEnergy()

    conf = mol_opt.GetConformer()
    coords = []
    for i in range(mol_opt.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])

    bonds = []
    for bond in mol_opt.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds.append({"i": i, "j": j, "order": bond.GetBondTypeAsDouble()})

    atoms = []
    for atom in mol_opt.GetAtoms():
        atoms.append(atom.GetSymbol())

    return {
        "name": name,
        "smiles": smiles,
        "n_atoms": mol.GetNumAtoms(),
        "atom_types": atom_types,
        "atoms": atoms,
        "final_energy": final_energy,
        "optimized_coords": coords,
        "bonds": bonds,
    }


def main():
    print("Gathering RDKit MMFF94s optimized geometries...")
    print("=" * 80)

    results = []
    for name, smiles in SMILES_LIST:
        print(f"  {name} ({smiles})...", end=" ", flush=True)
        r = get_rdkit_optimized(name, smiles)
        if r is None:
            print("SKIPPED")
            continue
        print(f"E_final={r['final_energy']:.4f} kcal/mol")
        results.append(r)

    json_path = "/tmp/rdkit_mmff94s_opt.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {json_path}")

    print("\n" + "=" * 80)
    print("SUMMARY: RDKit optimized geometries (these are the 'ground truth')")
    print("=" * 80)
    for r in results:
        print(f"\n  {r['name']}:")
        print(f"    Atoms: {r['atoms']}")
        print(f"    Atom types: {r['atom_types']}")
        print(f"    Energy: {r['final_energy']:.4f} kcal/mol")
        print(f"    Optimized coords:")
        for i, (c, a) in enumerate(zip(r["optimized_coords"], r["atoms"])):
            print(f"      {i} {a}: [{c[0]:.4f}, {c[1]:.4f}, {c[2]:.4f}]")
        print(f"    Bonds:")
        for b in r["bonds"]:
            i, j = b["i"], b["j"]
            p1, p2 = r["optimized_coords"][i], r["optimized_coords"][j]
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))
            print(f"      {i}({r['atoms'][i]})-{j}({r['atoms'][j]}): {d:.4f} A")


if __name__ == "__main__":
    main()
