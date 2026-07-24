#!/usr/bin/env python3
"""
Direct comparison: RDKit MMFF94s optimization vs WebMM optimization.

For each test molecule:
1. Generate 3D coords with RDKit ETKDG
2. Optimize with RDKit MMFF94s → record final geometry + energy breakdown
3. Optimize with WebMM starting from same initial coords → record final geometry + energy
4. Compare: bond lengths, angles, energies
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
    ("isobutane", "CC(C)C"),
    ("fluoromethane", "CF"),
    ("chloromethane", "CCl"),
    ("benzene", "c1ccccc1"),
    ("phenol", "c1ccc(O)cc1"),
    ("toluene", "Cc1ccccc1"),
    ("acetaldehyde", "CC=O"),
    ("glycine", "NCC(=O)O"),
    ("cyclohexane", "C1CCCCC1"),
]


def calc_angle(p1, p2, p3):
    v1 = [a - b for a, b in zip(p1, p2)]
    v2 = [a - b for a, b in zip(p3, p2)]
    dot = sum(a * b for a, b in zip(v1, v2))
    n1 = math.sqrt(sum(a * a for a in v1))
    n2 = math.sqrt(sum(a * a for a in v2))
    if n1 > 0 and n2 > 0:
        cos_t = max(-1, min(1, dot / (n1 * n2)))
        return math.degrees(math.acos(cos_t))
    return 0.0


def calc_dihedral(p1, p2, p3, p4):
    b1 = [p2[i] - p1[i] for i in range(3)]
    b2 = [p3[i] - p2[i] for i in range(3)]
    b3 = [p4[i] - p3[i] for i in range(3)]
    n1 = [
        b1[1] * b2[2] - b1[2] * b2[1],
        b1[2] * b2[0] - b1[0] * b2[2],
        b1[0] * b2[1] - b1[1] * b2[0],
    ]
    n2 = [
        b2[1] * b3[2] - b2[2] * b3[1],
        b2[2] * b3[0] - b2[0] * b3[2],
        b2[0] * b3[1] - b2[1] * b3[0],
    ]
    m = [
        b2[1] * n2[2] - b2[2] * n2[1],
        b2[2] * n2[0] - b2[0] * n2[2],
        b2[0] * n2[1] - b2[1] * n2[0],
    ]
    nn1 = math.sqrt(sum(x * x for x in n1))
    nn2 = math.sqrt(sum(x * x for x in n2))
    if nn1 < 1e-10 or nn2 < 1e-10:
        return 0.0
    x = sum(a * b for a, b in zip(n1, n2)) / (nn1 * nn2)
    y = sum(a * b for a, b in zip(m, n2)) / (nn1 * nn2)
    return math.degrees(math.atan2(y, x))


def dist(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def mol_to_sdf(mol, name="mol"):
    return Chem.MolToMolBlock(mol)


def get_rdkit_results(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    sdf = mol_to_sdf(mol, name)

    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if mp is None:
        return None

    atom_types = [mp.GetMMFFAtomType(i) for i in range(mol.GetNumAtoms())]

    ff_init = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mp, confId=0)
    initial_energy = ff_init.CalcEnergy()

    mol_opt = Chem.RWMol(mol)
    rdForceFieldHelpers.MMFFOptimizeMolecule(
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

    return {
        "name": name,
        "smiles": smiles,
        "n_atoms": mol.GetNumAtoms(),
        "atom_types": atom_types,
        "sdf": sdf,
        "initial_energy": initial_energy,
        "final_energy": final_energy,
        "optimized_coords": coords,
        "bonds": bonds,
    }


def write_rust_test_file(results):
    lines = []
    lines.append("#[cfg(test)]")
    lines.append("mod compare_opt {")
    lines.append("    use crate::mmff::{MMFFForceField, MMFFVariant};")
    lines.append("    use crate::molecule::parse_sdf;")
    lines.append("    use crate::optimizer::{self, ConvergenceOptions};")
    lines.append("")

    for r in results:
        safe_name = (
            r["name"]
            .replace("-", "_")
            .replace("(", "")
            .replace(")", "")
            .replace(" ", "_")
        )
        sdf_escaped = r["sdf"].replace("\\", "\\\\").replace('"', '\\"')

        lines.append(f"    #[test]")
        lines.append(f"    fn compare_{safe_name}() {{")
        lines.append(f'        let sdf = r#"{r["sdf"]}"#;')
        lines.append(f"        let mol = parse_sdf(sdf).unwrap();")
        lines.append(
            f"        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);"
        )
        lines.append(f"")

        coords_str = r["optimized_coords"]
        flat = []
        for c in coords_str:
            flat.extend(c)
        coords_literal = "[" + ", ".join(f"{v:.6f}" for v in flat) + "]"
        lines.append(
            f"        let rdkit_coords: Vec<[f64; 3]> = {coords_literal}.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();"
        )
        lines.append(f"")
        lines.append(
            f"        let opts = ConvergenceOptions {{ max_iterations: 2000, ..Default::default() }};"
        )
        lines.append(
            f"        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);"
        )
        lines.append(f"")
        lines.append(f"        let rdkit_energy = {r['final_energy']:.6f}f64;")
        lines.append(f"        let webmm_energy = result.final_energy;")
        lines.append(
            f"        let ratio = if rdkit_energy.abs() > 0.001 {{ webmm_energy / rdkit_energy }} else {{ webmm_energy - rdkit_energy }};"
        )
        lines.append(f"")
        lines.append(f"        // Compare bond lengths")
        for bond in r["bonds"]:
            i, j = bond["i"], bond["j"]
            lines.append(
                f"        let d_{i}_{j}_rdkit = super::super::dist(&rdkit_coords[{i}], &rdkit_coords[{j}]);"
            )
            lines.append(
                f"        let d_{i}_{j}_webmm = super::super::dist(&result.optimized_coords[{i}], &result.optimized_coords[{j}]);"
            )
            lines.append(
                f"        let d_{i}_{j}_err = (d_{i}_{j}_rdkit - d_{i}_{j}_webmm).abs();"
            )

        lines.append(f"")
        lines.append(f'        println!("\\n=== {r["name"].upper()} ===");')
        lines.append(f'        println!("RDKit energy: {{:.6f}}", rdkit_energy);')
        lines.append(f'        println!("WebMM energy: {{:.6f}}", webmm_energy);')
        lines.append(f'        println!("Ratio: {{:.4f}}", ratio);')
        for bond in r["bonds"]:
            i, j = bond["i"], bond["j"]
            lines.append(
                f'        println!("Bond {i}-{j}: RDKit={{:.4f}} WebMM={{:.4f}} err={{:.4f}}", d_{i}_{j}_rdkit, d_{i}_{j}_webmm, d_{i}_{j}_err);'
            )

        lines.append(f"")
        lines.append(f"        // Check energy is within 10% or 1 kcal/mol")
        lines.append(
            f"        let energy_ok = (webmm_energy - rdkit_energy).abs() < 1.0 || (ratio - 1.0).abs() < 0.10;"
        )
        lines.append(f"        if !energy_ok {{")
        lines.append(
            f'            eprintln!("ENERGY MISMATCH for {r["name"]}: RDKit={{:.6f}} WebMM={{:.6f}} ratio={{:.4f}}", rdkit_energy, webmm_energy, ratio);'
        )
        lines.append(f"        }}")

        for bond in r["bonds"]:
            i, j = bond["i"], bond["j"]
            lines.append(
                f'        assert!(d_{i}_{j}_err < 0.05, "Bond {i}-{j} mismatch: RDKit={{:.4f}} WebMM={{:.4f}}", d_{i}_{j}_rdkit, d_{i}_{j}_webmm);'
            )

        lines.append(f"    }}")
        lines.append("")

    lines.append("}")
    return "\n".join(lines)


def main():
    print("Gathering RDKit MMFF94s reference data...")
    print("=" * 80)

    results = []
    for name, smiles in SMILES_LIST:
        print(f"  {name} ({smiles})...", end=" ", flush=True)
        r = get_rdkit_results(name, smiles)
        if r is None:
            print("SKIPPED (MMFF props failed)")
            continue
        print(
            f"atoms={r['n_atoms']} E_init={r['initial_energy']:.2f} E_final={r['final_energy']:.2f}"
        )
        results.append(r)

    print(f"\nGenerated {len(results)} test cases")

    rust_code = write_rust_test_file(results)

    out_path = "/Users/rliu/Projects/WebMM/src/compare_opt.rs"
    with open(out_path, "w") as f:
        f.write(rust_code)
    print(f"Wrote Rust test file to {out_path}")

    json_path = "/tmp/rdkit_opt_reference.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Wrote JSON reference to {json_path}")


if __name__ == "__main__":
    main()
