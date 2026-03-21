//! SDF/MOL file parser

use super::{Atom, Bond, BondType, Molecule};

/// Parse SDF/MOL file content into Molecule structure
pub fn parse_sdf(sdf_content: &str) -> Result<Molecule, String> {
    let lines: Vec<&str> = sdf_content.lines().collect();

    if lines.len() < 4 {
        return Err("SDF file too short".to_string());
    }

    // Find counts line (typically line 3-4, but search for line with counts)
    let counts_line_idx = lines
        .iter()
        .position(|line| {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.len() < 10 {
                return false;
            }
            let first_10 = trimmed.chars().take(10).collect::<String>();
            let parts: Vec<&str> = first_10.split_whitespace().collect();
            if parts.len() < 3 {
                return false;
            }
            match parts[0].parse::<usize>() {
                Ok(count) => count > 0 && count <= 1000,
                Err(_) => false,
            }
        })
        .ok_or_else(|| "Could not find counts line".to_string())?;

    let counts_line = lines[counts_line_idx].trim();
    let parts: Vec<&str> = counts_line.split_whitespace().collect();

    if parts.len() < 3 {
        return Err("Invalid counts line".to_string());
    }

    let num_atoms: usize = parts[0]
        .trim()
        .parse()
        .map_err(|e| format!("Failed to parse atom count: {}", e))?;

    let num_bonds: usize = parts[1]
        .trim()
        .parse()
        .map_err(|e| format!("Failed to parse bond count: {}", e))?;

    // Parse atom block (after counts line)
    let atom_start = counts_line_idx + 1;
    let atom_end = atom_start + num_atoms;

    if lines.len() < atom_end {
        return Err(format!(
            "Not enough atom lines: expected {}, got {}",
            atom_end,
            lines.len()
        ));
    }

    let mut atoms = Vec::with_capacity(num_atoms);
    let mut bonds = Vec::new();

    for (i, line) in lines[atom_start..atom_end].iter().enumerate() {
        if line.len() < 34 {
            return Err(format!("Atom line {} too short: {}", i + 1, line.len()));
        }

        let end = std::cmp::min(34, line.len());
        let symbol = line[31..end].trim();
        let symbol = if symbol.is_empty() {
            "C".to_string()
        } else {
            symbol.to_string()
        };

        // Extract coordinates (columns 0-9, 10-19, 20-29)
        let x: f64 = line[0..9].trim().parse().unwrap_or(0.0);
        let y: f64 = line[10..19].trim().parse().unwrap_or(0.0);
        let z: f64 = line[20..29].trim().parse().unwrap_or(0.0);

        let atomic_number = get_atomic_number(&symbol);
        let mass = get_atomic_mass(atomic_number);
        let charge = parse_charge(line);

        atoms.push(Atom {
            symbol,
            atomic_number,
            mass,
            charge,
            position: [x, y, z],
            index: i,
        });
    }

    // Parse bond block (after atom block)
    let bond_start = atom_end;
    let bond_end = bond_start + num_bonds;

    if lines.len() >= bond_end {
        for line in lines[bond_start..bond_end].iter() {
            if line.len() < 7 {
                continue;
            }

            let atom1: usize = line[0..3]
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("Failed to parse atom1: {}", e))?
                .checked_sub(1)
                .ok_or_else(|| "Invalid atom1 index".to_string())?;

            let atom2: usize = line[3..6]
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("Failed to parse atom2: {}", e))?
                .checked_sub(1)
                .ok_or_else(|| "Invalid atom2 index".to_string())?;

            let bond_type_num: i32 = line[6..9]
                .trim()
                .parse::<i32>()
                .map_err(|e| format!("Failed to parse bond type: {}", e))?;

            let bond_type = match bond_type_num {
                1 => BondType::Single,
                2 => BondType::Double,
                3 => BondType::Triple,
                4 => BondType::Aromatic,
                5..=7 => BondType::Aromatic,
                _ => BondType::Single,
            };

            bonds.push(Bond {
                atom1,
                atom2,
                bond_type,
            });
        }
    }

    // Extract molecule name from header
    // V2000 format: line 0 = name, line 1 = program/timestamp, line 2 = comment, line 3 = counts
    let name = if counts_line_idx >= 3 {
        lines[counts_line_idx - 3].trim().to_string()
    } else {
        lines[0].trim().to_string()
    };

    let adjacency = super::graph::build_adjacency_list_from_bonds(&atoms.len(), &bonds);

    Ok(Molecule {
        atoms,
        bonds,
        name,
        adjacency,
    })
}

/// Get atomic number from element symbol
fn get_atomic_number(symbol: &str) -> u8 {
    match symbol.trim() {
        "H" => 1,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Br" => 35,
        "I" => 53,
        _ => 0,
    }
}

/// Get atomic mass in atomic mass units (amu)
fn get_atomic_mass(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 1.00794,
        6 => 12.0107,
        7 => 14.0067,
        8 => 15.9994,
        9 => 18.9984,
        15 => 30.9738,
        16 => 32.0660,
        17 => 35.4527,
        35 => 79.9040,
        53 => 126.9045,
        _ => 12.0107,
    }
}

/// Parse charge from atom line using V2000 encoding
fn parse_charge(line: &str) -> f64 {
    // V2000 charge field encoding (column 36-39)
    // 0=none, 1=+3, 2=+2, 3=+1, 4=doublet radical, 5=-1, 6=-2, 7=-3
    if line.len() > 39 {
        let charge_code: i32 = line[36..40].trim().parse().unwrap_or(0);
        match charge_code {
            0 => 0.0,
            1 => 3.0,
            2 => 2.0,
            3 => 1.0,
            4 => 0.0,
            5 => -1.0,
            6 => -2.0,
            7 => -3.0,
            _ => 0.0,
        }
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_water_sdf() {
        let sdf = "Water
     CDK     101218203532D 0

   3  2  0  0  0  0            999 V2000
     0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
     0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END";

        let mol = parse_sdf(sdf).expect("Failed to parse water SDF");
        assert_eq!(mol.name, "Water");
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.bonds.len(), 2);

        assert_eq!(mol.atoms[0].symbol, "O");
        assert_eq!(mol.atoms[0].atomic_number, 8);
        assert_eq!(mol.atoms[1].symbol, "H");
        assert_eq!(mol.atoms[1].atomic_number, 1);
        assert_eq!(mol.atoms[2].symbol, "H");
        assert_eq!(mol.atoms[2].atomic_number, 1);

        assert_eq!(mol.adjacency.len(), 3);
        assert_eq!(mol.adjacency[0], vec![1, 2]);
        assert_eq!(mol.adjacency[1], vec![0]);
        assert_eq!(mol.adjacency[2], vec![0]);
    }

    #[test]
    fn test_parse_charge_encoding() {
        let sdf = "Charged
     CDK     101218203532D 0

   1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O  0  5  0  0  0  0  0  0  0  0  0  0  0
M  END";

        let mol = parse_sdf(sdf).expect("Failed to parse charged SDF");
        assert_eq!(mol.atoms[0].charge, -1.0);
    }
}
