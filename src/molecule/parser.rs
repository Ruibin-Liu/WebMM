//! SDF/MOL file parser

use super::{Atom, Bond, BondStereo, BondType, Molecule};

/// Parse SDF/MOL file content into Molecule structure
/// Supports both V2000 and V3000 formats
pub fn parse_sdf(sdf_content: &str) -> Result<Molecule, String> {
    let lines: Vec<&str> = sdf_content.lines().collect();

    if lines.len() < 4 {
        return Err("SDF file too short".to_string());
    }

    // Detect format version
    let is_v3000 = lines.iter().any(|line| line.contains("V3000"));

    if is_v3000 {
        parse_v3000(&lines)
    } else {
        parse_v2000(&lines)
    }
}

/// Parse V2000 format MOL/SDF file
fn parse_v2000(lines: &[&str]) -> Result<Molecule, String> {
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
        // Extract coordinates using whitespace splitting (more robust than fixed columns)
        let prefix = line.get(..40).unwrap_or(line);
        let parts: Vec<&str> = prefix.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(format!("Atom line {} has too few fields: {}", i + 1, line));
        }
        let x: f64 = parts[0].parse().unwrap_or(0.0);
        let y: f64 = parts[1].parse().unwrap_or(0.0);
        let z: f64 = parts[2].parse().unwrap_or(0.0);
        let symbol = parts[3].to_string();

        let atomic_number = get_atomic_number(&symbol);
        let mass = get_atomic_mass(atomic_number);
        let charge = parse_charge(line);
        let stereo_parity = parse_stereo_parity(line);

        atoms.push(Atom {
            symbol,
            atomic_number,
            mass,
            charge,
            position: [x, y, z],
            index: i,
            stereo_parity,
        });
    }

    // Parse bond block (after atom block)
    let bond_start = atom_end;
    let bond_end = bond_start + num_bonds;

    if lines.len() >= bond_end {
        for line in lines[bond_start..bond_end].iter() {
            // A well-formed V2000 bond line carries three 3-column fields
            // (18 chars with columns). Shorter lines are skipped (historic
            // lenient behavior); 7-8 char lines used to PANIC on the
            // [6..9] slice below (fixed-column get is boundary-safe).
            if line.len() < 9 {
                continue;
            }

            // get(a..b) is None on out-of-bounds OR a non-char-boundary
            // (multi-byte UTF-8) — both are malformed fixed-column fields
            let f1 = line
                .get(0..3)
                .ok_or_else(|| format!("Malformed bond line (atom1 field): {:?}", line))?;
            let atom1: usize = f1
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("Failed to parse atom1: {}", e))?
                .checked_sub(1)
                .ok_or_else(|| "Invalid atom1 index".to_string())?;

            let f2 = line
                .get(3..6)
                .ok_or_else(|| format!("Malformed bond line (atom2 field): {:?}", line))?;
            let atom2: usize = f2
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("Failed to parse atom2: {}", e))?
                .checked_sub(1)
                .ok_or_else(|| "Invalid atom2 index".to_string())?;

            // Out-of-range references must be REJECTED, not silently
            // dropped: mol.bonds would disagree with atoms/adjacency and
            // every bonds-iterating consumer (ring perception, typing)
            // indexes per-atom arrays with these (was a panic chain)
            if atom1 >= num_atoms || atom2 >= num_atoms {
                return Err(format!(
                    "Bond {}-{} references an atom out of range (molecule has {} atoms)",
                    atom1 + 1,
                    atom2 + 1,
                    num_atoms
                ));
            }

            let f3 = line
                .get(6..9)
                .ok_or_else(|| format!("Malformed bond line (type field): {:?}", line))?;
            let bond_type_num: i32 = f3
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

            let stereo = if line.len() >= 12 {
                // boundary-safe: None on short or multi-byte-column lines
                let stereo_code: i32 = line
                    .get(9..12)
                    .and_then(|s| s.trim().parse().ok())
                    .unwrap_or(0);
                match stereo_code {
                    // V2000: 1 = wedge up, 6 = hash down (CTFile spec)
                    1 => BondStereo::Wedge,
                    6 => BondStereo::Hash,
                    3 => BondStereo::Cis,
                    7 => BondStereo::Trans,
                    _ => BondStereo::None,
                }
            } else {
                BondStereo::None
            };

            bonds.push(Bond {
                atom1,
                atom2,
                bond_type,
                stereo,
                kekule_type: None,
            });
        }
    }

    // Parse M  CHG property records (authoritative formal charges in modern
    // V2000 SDFs, e.g. all RDKit-generated files; atom-line charge codes are 0)
    if lines.len() >= bond_end {
        for line in lines[bond_end..].iter() {
            let trimmed = line.trim_end();
            if trimmed.starts_with("M  END") {
                break;
            }
            if let Some(rest) = trimmed.strip_prefix("M  CHG") {
                let nums: Vec<i64> = rest
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                let Some(&count) = nums.first() else { continue };
                for pair in nums[1..].chunks(2).take(count as usize) {
                    if let [atom_num, chg] = pair {
                        if *atom_num >= 1 && (*atom_num as usize) <= atoms.len() {
                            atoms[(*atom_num - 1) as usize].charge = *chg as f64;
                        }
                    }
                }
            }
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

    let mut mol = Molecule {
        atoms,
        bonds,
        name,
        adjacency,
    };
    super::graph::perceive_aromatic_bonds(&mut mol);
    Ok(mol)
}

/// Parse V3000 format MOL/SDF file
fn parse_v3000(lines: &[&str]) -> Result<Molecule, String> {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut name = String::new();
    let mut in_atom_block = false;
    let mut in_bond_block = false;

    // First line after header is molecule name
    if !lines.is_empty() {
        name = lines[0].trim().to_string();
    }

    for line in lines.iter() {
        let trimmed = line.trim();

        // Handle BEGIN/END blocks (be flexible with spacing - V3000 uses "M  V30")
        if trimmed.contains("BEGIN ATOM") && trimmed.starts_with("M") && trimmed.contains("V30") {
            in_atom_block = true;
            continue;
        }
        if trimmed.contains("END ATOM") && trimmed.starts_with("M") && trimmed.contains("V30") {
            in_atom_block = false;
            continue;
        }
        if trimmed.contains("BEGIN BOND") && trimmed.starts_with("M") && trimmed.contains("V30") {
            in_bond_block = true;
            continue;
        }
        if trimmed.contains("END BOND") && trimmed.starts_with("M") && trimmed.contains("V30") {
            in_bond_block = false;
            continue;
        }

        // Skip non-data lines (including BEGIN/END CTAB, COUNTS, etc.)
        if !(trimmed.starts_with("M") && trimmed.contains("V30")) {
            continue;
        }

        // Skip non-atom/bond lines
        if trimmed.contains("BEGIN") || trimmed.contains("END") || trimmed.contains("COUNTS") {
            continue;
        }

        // Remove the "M  V30" prefix (V3000 format uses two spaces after M)
        let data = if let Some(pos) = trimmed.find("V30") {
            trimmed[pos + 3..].trim()
        } else {
            trimmed
        };

        if in_atom_block {
            // Atom format: index symbol x y z [optional fields]
            let parts: Vec<&str> = data.split_whitespace().collect();
            if parts.len() >= 5 {
                let index: usize = parts[0].parse().map_err(|_| "Invalid atom index")?;
                let symbol = parts[1].to_string();
                let x: f64 = parts[2].parse().map_err(|_| "Invalid x coordinate")?;
                let y: f64 = parts[3].parse().map_err(|_| "Invalid y coordinate")?;
                let z: f64 = parts[4].parse().map_err(|_| "Invalid z coordinate")?;

                let atomic_number = get_atomic_number(&symbol);
                let mass = get_atomic_mass(atomic_number);
                let charge = 0.0; // V3000 stores charge separately
                let stereo_parity = parse_v3000_stereo_parity(data);

                atoms.push(Atom {
                    symbol,
                    atomic_number,
                    mass,
                    charge,
                    position: [x, y, z],
                    index: index - 1, // Convert to 0-based
                    stereo_parity,
                });
            }
        } else if in_bond_block {
            // Bond format: index type atom1 atom2
            let parts: Vec<&str> = data.split_whitespace().collect();
            if parts.len() >= 4 {
                let bond_type_num: i32 = parts[1].parse().map_err(|_| "Invalid bond type")?;
                let atom1: usize = parts[2].parse::<usize>().map_err(|_| "Invalid atom1")?;
                let atom2: usize = parts[3].parse::<usize>().map_err(|_| "Invalid atom2")?;

                let bond_type = match bond_type_num {
                    1 => BondType::Single,
                    2 => BondType::Double,
                    3 => BondType::Triple,
                    4 => BondType::Aromatic,
                    _ => BondType::Single,
                };

                bonds.push(Bond {
                    atom1: atom1 - 1, // Convert to 0-based
                    atom2: atom2 - 1,
                    bond_type,
                    ..Default::default()
                });
            }
        }
    }

    // Sort atoms by index to ensure correct order
    atoms.sort_by_key(|a| a.index);

    // Re-index atoms to be 0-based sequential
    for (i, atom) in atoms.iter_mut().enumerate() {
        atom.index = i;
    }

    let adjacency = super::graph::build_adjacency_list_from_bonds(&atoms.len(), &bonds);

    let mut mol = Molecule {
        atoms,
        bonds,
        name,
        adjacency,
    };
    super::graph::perceive_aromatic_bonds(&mut mol);
    Ok(mol)
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
        "Li" => 3,
        "Na" => 11,
        "K" => 19,
        "Mg" => 12,
        "Ca" => 20,
        "Zn" => 30,
        "Fe" => 26,
        "Cu" => 29,
        "Si" => 14,
        // full first rows for metal complexes (GFN-FF covers H-Zn + Br/I)
        "He" => 2,
        "Be" => 4,
        "B" => 5,
        "Ne" => 10,
        "Al" => 13,
        "Ar" => 18,
        "Sc" => 21,
        "Ti" => 22,
        "V" => 23,
        "Cr" => 24,
        "Mn" => 25,
        "Co" => 27,
        "Ni" => 28,
        "Ga" => 31,
        "Ge" => 32,
        "As" => 33,
        "Se" => 34,
        "Kr" => 36,
        "Sn" => 50,
        "Sb" => 51,
        "Te" => 52,
        "Xe" => 54,
        "Pb" => 82,
        "Bi" => 83,
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
        3 => 6.941,
        11 => 22.9898,
        19 => 39.0983,
        12 => 24.305,
        20 => 40.078,
        30 => 65.38,
        26 => 55.845,
        29 => 63.546,
        14 => 28.0855,
        _ => 12.0107,
    }
}

/// Parse charge from atom line using V2000 encoding
fn parse_charge(line: &str) -> f64 {
    // V2000 charge field encoding (column 36-39)
    // 0=none, 1=+3, 2=+2, 3=+1, 4=doublet radical, 5=-1, 6=-2, 7=-3
    if line.len() > 39 {
        // boundary-safe on non-ASCII lines (get is None → charge 0, was a panic)
        let charge_code: i32 = line
            .get(36..40)
            .map(|f| f.trim().parse().unwrap_or(0))
            .unwrap_or(0);
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

fn parse_stereo_parity(line: &str) -> u8 {
    // V2000 atom stereo parity field (column 39-42)
    // 0=not specified, 1=odd (CCW), 2=even (CW), 3=either
    if line.len() > 42 {
        line.get(39..42)
            .map(|f| f.trim().parse().unwrap_or(0))
            .unwrap_or(0)
    } else {
        0
    }
}

fn parse_v3000_stereo_parity(data: &str) -> u8 {
    // V3000: CFG=n field in atom line
    // 1=odd (CCW/R), 2=even (CW/S), 3=either
    for part in data.split_whitespace() {
        if let Some(val) = part.strip_prefix("CFG=") {
            let cfg: u8 = val.parse().unwrap_or(0);
            return match cfg {
                1 => 1,
                2 => 2,
                3 => 3,
                _ => 0,
            };
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- review regressions: malformed input must Err/skip, never panic ----
    const PROBE_HDR: &str = "probe\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n";
    const C_LINE: &str =
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n";
    const C_LINE2: &str =
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n";

    #[test]
    fn bond_line_short_is_skipped_not_panicking() {
        // 8-char bond line: used to panic on the [6..9] slice
        let sdf = format!("{}{}{} 1  2  1\nM  END\n", PROBE_HDR, C_LINE, C_LINE2);
        let r = parse_sdf(&sdf);
        assert!(
            r.is_ok(),
            "short bond line should be skipped, got {:?}",
            r.err()
        );
        assert_eq!(r.unwrap().bonds.len(), 0);
    }

    #[test]
    fn bond_out_of_range_is_rejected() {
        // used to be silently kept in mol.bonds → inconsistent state →
        // downstream per-atom indexing panics (ring perception, typing)
        let sdf = format!(
            "{}{}{}  1  9  1  0  0  0\nM  END\n",
            PROBE_HDR, C_LINE, C_LINE2
        );
        let err = parse_sdf(&sdf).unwrap_err();
        assert!(err.contains("out of range"), "got: {}", err);
    }

    #[test]
    fn multibyte_columns_do_not_panic() {
        // non-ASCII crossing a fixed-column boundary used to panic the slice
        let atom_line = "    0.0000    0.0000    0.0000 C  µ²0  0  0  0  0  0  0  0  0\n";
        let sdf = format!(
            "{}{}{}  1  2  1  0  0  0\nM  END\n",
            PROBE_HDR, atom_line, C_LINE2
        );
        let _ = parse_sdf(&sdf); // must not panic
        let bond_line = "  1 µ2  1  0  0  0\n";
        let sdf = format!("{}{}{}{}M  END\n", PROBE_HDR, C_LINE, C_LINE2, bond_line);
        let _ = parse_sdf(&sdf); // must not panic (Err is fine)
    }

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

    #[test]
    fn test_parse_v3000_water() {
        let sdf = "Water
  CDK

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.0000 0.0000 0.0000 0
M  V30 2 H 0.9580 0.0000 0.0000 0
M  V30 3 H -0.2390 0.9270 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 END BOND
M  V30 END CTAB
M  END";

        let mol = parse_sdf(sdf).expect("Failed to parse V3000 water SDF");
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
    fn test_parse_v3000_methane() {
        let sdf = "Methane
  CDK

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 H 0.6290 0.6290 0.6290 0
M  V30 3 H 0.6290 -0.6290 -0.6290 0
M  V30 4 H -0.6290 0.6290 -0.6290 0
M  V30 5 H -0.6290 -0.6290 0.6290 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END";

        let mol = parse_sdf(sdf).expect("Failed to parse V3000 methane SDF");
        assert_eq!(mol.name, "Methane");
        assert_eq!(mol.atoms.len(), 5);
        assert_eq!(mol.bonds.len(), 4);

        assert_eq!(mol.atoms[0].symbol, "C");
        assert_eq!(mol.atoms[0].atomic_number, 6);
    }
}
