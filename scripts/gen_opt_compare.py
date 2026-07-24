#!/usr/bin/env python3
"""Generate opt_compare.rs with RDKit reference data."""

import json

with open("/tmp/rdkit_mmff94s_opt.json") as f:
    data = json.load(f)

lines = []
lines.append("#[cfg(test)]")
lines.append("mod opt_compare {")
lines.append("    use crate::mmff::{MMFFForceField, MMFFVariant};")
lines.append("    use crate::molecule::{Molecule, Atom, Bond, BondType, BondStereo};")
lines.append("    use crate::optimizer;")
lines.append("    use crate::ConvergenceOptions;")
lines.append("")
lines.append("    fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {")
lines.append(
    "        ((a[0]-b[0]).powi(2)+(a[1]-b[1]).powi(2)+(a[2]-b[2]).powi(2)).sqrt()"
)
lines.append("    }")
lines.append("")
lines.append(
    "    fn make_mol(name: &str, atom_syms: &[&str], bonds: &[(usize, usize, u8)], coords: &[[f64; 3]]) -> Molecule {"
)
lines.append(
    '        let anum = |s: &str| -> u8 { match s { "H"=>1,"C"=>6,"N"=>7,"O"=>8,"F"=>9,"P"=>15,"S"=>16,"Cl"=>17,"Br"=>35,"I"=>53,_=>0 } };'
)
lines.append(
    "        let atoms: Vec<Atom> = atom_syms.iter().enumerate().map(|(i, s)| Atom {"
)
lines.append("            symbol: s.to_string(),")
lines.append("            atomic_number: anum(s),")
lines.append("            mass: 0.0,")
lines.append("            charge: 0.0,")
lines.append("            position: coords[i],")
lines.append("            index: i,")
lines.append("            stereo_parity: 0,")
lines.append("        }).collect();")
lines.append("        let bonds: Vec<Bond> = bonds.iter().map(|(a, b, bo)| Bond {")
lines.append("            atom1: *a, atom2: *b,")
lines.append(
    "            bond_type: match bo { 20 => BondType::Double, 30 => BondType::Triple, 15 => BondType::Aromatic, _ => BondType::Single },"
)
lines.append("            stereo: BondStereo::None,")
lines.append("        }).collect();")
lines.append("        let n = atoms.len();")
lines.append("        let mut adj = vec![vec![]; n];")
lines.append(
    "        for b in &bonds { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }"
)
lines.append(
    "        Molecule { atoms, bonds, name: name.to_string(), adjacency: adj }"
)
lines.append("    }")

for mol in data:
    name = (
        mol["name"]
        .replace("-", "_")
        .replace("(", "")
        .replace(")", "")
        .replace(" ", "_")
    )
    rdkit_e = mol["final_energy"]
    coords = mol["optimized_coords"]
    bonds = mol["bonds"]
    atoms = mol["atoms"]
    display = mol["name"]

    syms = ", ".join(f'"{a}"' for a in atoms)
    bond_strs = ", ".join(
        f"({b['i']}, {b['j']}, {int(b['order'] * 10)})" for b in bonds
    )
    coord_strs = ", ".join(f"[{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]" for c in coords)

    lines.append("")
    lines.append("    #[test]")
    lines.append(f"    fn opt_{name}() {{")
    lines.append(f"        let syms: &[&str] = &[{syms}];")
    lines.append(f"        let bonds: &[(usize, usize, u8)] = &[{bond_strs}];")
    lines.append(f"        let coords: &[[f64; 3]] = &[{coord_strs}];")
    lines.append(f'        let mol = make_mol("{display}", syms, bonds, coords);')
    lines.append("        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);")
    lines.append("        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();")
    lines.append(
        "        let opts = ConvergenceOptions { max_iterations: 2000, ..Default::default() };"
    )
    lines.append("        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);")
    lines.append(f"        let rdkit_e = {rdkit_e:.6f}_f64;")
    lines.append("        let webmm_e = result.final_energy;")
    lines.append(f'        println!("\\n=== {display.upper()} ===");')
    lines.append(
        '        println!("RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}", rdkit_e, webmm_e, webmm_e - rdkit_e);'
    )
    lines.append(
        '        println!("Converged: {} ({} iters)", result.converged, result.iterations);'
    )

    for b in bonds:
        i, j = b["i"], b["j"]
        lines.append("        {")
        lines.append(
            f"            let dr = dist(&rdkit_coords[{i}], &rdkit_coords[{j}]);"
        )
        lines.append(
            f"            let dw = dist(&result.optimized_coords[{i}], &result.optimized_coords[{j}]);"
        )
        lines.append("            let err = (dr - dw).abs();")
        lines.append(
            f'            if err > 0.01 {{ println!("  Bond {i}-{j}: RDKit={{:.4}} WebMM={{:.4}} err={{:.4}}", dr, dw, err); }}'
        )
        lines.append("        }")

    lines.append("    }")

lines.append("}")

with open("/Users/rliu/Projects/WebMM/src/opt_compare.rs", "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"Generated {len(data)} tests")
