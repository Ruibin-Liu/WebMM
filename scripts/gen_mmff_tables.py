"""Generate corrected mmff_tables.rs with proper EQ levels, angle, and SB tables."""

import sys


def parse_def():
    entries = {}
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFDef.tsv"
    ) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                try:
                    at = int(parts[1])
                    entries[at] = [
                        int(parts[2]),
                        int(parts[3]),
                        int(parts[4]),
                        int(parts[5]),
                    ]
                except:
                    pass
    return entries


def parse_angle():
    entries = []
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFAngle.tsv"
    ) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                try:
                    entries.append(
                        (
                            int(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            int(parts[3]),
                            float(parts[4]),
                            float(parts[5]),
                        )
                    )
                except:
                    pass
    return entries


def parse_stbn():
    entries = []
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFStbn.tsv"
    ) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                try:
                    entries.append(
                        (
                            int(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            int(parts[3]),
                            float(parts[4]),
                            float(parts[5]),
                        )
                    )
                except:
                    pass
    return entries


def parse_dfsb():
    entries = []
    with open(
        "/Users/rliu/Projects/WebMM/mmff_params_extracted/defaultMMFFDfsb.tsv"
    ) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                try:
                    entries.append(
                        (
                            int(parts[0]),
                            int(parts[1]),
                            int(parts[2]),
                            float(parts[3]),
                            float(parts[4]),
                        )
                    )
                except:
                    pass
    return entries


def fmt_f64(v):
    s = f"{v:.6}"
    # Remove trailing zeros but keep at least one decimal
    if "." in s:
        s = s.rstrip("0")
        if s.endswith("."):
            s += "0"
    return s


out = []

# EQ levels table
def_data = parse_def()
out.append("#![allow(clippy::too_many_arguments)]")
out.append("")
out.append("#[rustfmt::skip]")
out.append("pub const MMFF_DEF_EQ_LEVELS: [[u8; 4]; 100] = [")
for i in range(100):
    if i in def_data:
        l = def_data[i]
        out.append(f"    [{l[0]:2}, {l[1]:2}, {l[2]:2}, {l[3]:2}],")
    else:
        out.append(f"    [ 0,  0,  0,  0],")
out.append("];")
out.append("")

# Angle table
angle_data = parse_angle()
out.append(f"#[rustfmt::skip]")
out.append(
    f"pub const MMFF_ANGLE_TABLE: [(u8, u8, u8, u8, f64, f64); {len(angle_data)}] = ["
)
for at, i, j, k, ka, t0 in angle_data:
    out.append(f"    ({at}, {i}, {j}, {k}, {fmt_f64(ka)}, {fmt_f64(t0)}),")
out.append("];")
out.append("")

# Stretch-bend table
stbn_data = parse_stbn()
out.append(f"#[rustfmt::skip]")
out.append(
    f"pub const MMFF_STBN_TABLE: [(u8, u8, u8, u8, f64, f64); {len(stbn_data)}] = ["
)
for sbt, i, j, k, kba_ijk, kba_kji in stbn_data:
    out.append(f"    ({sbt}, {i}, {j}, {k}, {fmt_f64(kba_ijk)}, {fmt_f64(kba_kji)}),")
out.append("];")
out.append("")

# Default stretch-bend table
dfsb_data = parse_dfsb()
out.append(f"#[rustfmt::skip]")
out.append(f"pub const MMFF_DFSB_TABLE: [(u8, u8, u8, f64, f64); {len(dfsb_data)}] = [")
for ir, jr, kr, kba_ijk, kba_kji in dfsb_data:
    out.append(f"    ({ir}, {jr}, {kr}, {fmt_f64(kba_ijk)}, {fmt_f64(kba_kji)}),")
out.append("];")
out.append("")

# Periodic table row lookup
out.append("pub fn get_periodic_table_row(atomic_number: u8) -> u8 {")
out.append("    match atomic_number {")
out.append("        1 => 0,  // H - period 1")
out.append("        6..=9 => 1,  // C, N, O, F - period 2")
out.append("        14..=17 => 2,  // Si, P, S, Cl - period 3")
out.append("        35 => 3,  // Br - period 4")
out.append("        53 => 4,  // I - period 5")
out.append("        _ => 1,")
out.append("    }")
out.append("}")
out.append("")

# Angle type computation
out.append(
    "pub fn compute_angle_type(bond_type_ij: u8, bond_type_jk: u8, ring_size: u8) -> u8 {"
)
out.append("    let bond_type_sum = bond_type_ij + bond_type_jk;")
out.append("    let mut angle_type = bond_type_sum;")
out.append("    if ring_size == 3 || ring_size == 4 {")
out.append("        angle_type = ring_size;")
out.append("        if bond_type_sum > 0 {")
out.append("            angle_type += bond_type_sum + ring_size - 2;")
out.append("        }")
out.append("    }")
out.append("    angle_type")
out.append("}")
out.append("")

# Stretch-bend type computation
out.append(
    "pub fn compute_stretch_bend_type(angle_type: u8, bond_type_1: u8, bond_type_2: u8) -> u8 {"
)
out.append("    match angle_type {")
out.append("        1 => {")
out.append(
    "            if bond_type_1 != 0 || bond_type_1 == bond_type_2 { 1 } else { 2 }"
)
out.append("        }")
out.append("        2 => 3,")
out.append("        4 => 4,")
out.append("        3 => 5,")
out.append("        5 => {")
out.append(
    "            if bond_type_1 != 0 || bond_type_1 == bond_type_2 { 6 } else { 7 }"
)
out.append("        }")
out.append("        6 => 8,")
out.append("        7 => {")
out.append(
    "            if bond_type_1 != 0 || bond_type_1 == bond_type_2 { 9 } else { 10 }"
)
out.append("        }")
out.append("        8 => 11,")
out.append("        _ => 0,")
out.append("    }")
out.append("}")
out.append("")

# Multi-stage angle lookup
out.append("pub fn lookup_angle_params(")
out.append("    type_i: u8,")
out.append("    type_j: u8,")
out.append("    type_k: u8,")
out.append("    angle_type: u8,")
out.append(") -> Option<(f64, f64)> {")
out.append("    let mut iter = 0;")
out.append("    while iter < 4 {")
out.append("        let idx_i = type_i as usize;")
out.append("        let idx_k = type_k as usize;")
out.append(
    "        let mut can_i = MMFF_DEF_EQ_LEVELS.get(idx_i).map_or(0, |l| l[iter]);"
)
out.append(
    "        let mut can_k = MMFF_DEF_EQ_LEVELS.get(idx_k).map_or(0, |l| l[iter]);"
)
out.append("        if can_i > can_k {")
out.append("            std::mem::swap(&mut can_i, &mut can_k);")
out.append("        }")
out.append("        if can_i == 0 && can_k == 0 && iter > 0 {")
out.append("            break;")
out.append("        }")
out.append("        for &(at, ti, tj, tk, ka, t0) in MMFF_ANGLE_TABLE.iter() {")
out.append(
    "            if at == angle_type && ti == can_i && tj == type_j && tk == can_k {"
)
out.append("                return Some((ka, t0));")
out.append("            }")
out.append("        }")
out.append("        iter += 1;")
out.append("    }")
out.append("    None")
out.append("}")
out.append("")

# Stretch-bend lookup
out.append("pub fn lookup_stretch_bend_params(")
out.append("    type_i: u8,")
out.append("    type_j: u8,")
out.append("    type_k: u8,")
out.append("    sb_type: u8,")
out.append("    atomic_num_i: u8,")
out.append("    atomic_num_j: u8,")
out.append("    atomic_num_k: u8,")
out.append(") -> Option<(f64, f64)> {")
out.append("    let mut can_i = type_i;")
out.append("    let mut can_k = type_k;")
out.append("    let mut swap = false;")
out.append("    if can_i > can_k {")
out.append("        std::mem::swap(&mut can_i, &mut can_k);")
out.append("        swap = true;")
out.append("    }")
out.append("    for &(st, ti, tj, tk, kba_ijk, kba_kji) in MMFF_STBN_TABLE.iter() {")
out.append("        if st == sb_type && ti == can_i && tj == type_j && tk == can_k {")
out.append("            if swap {")
out.append("                return Some((kba_kji, kba_ijk));")
out.append("            } else {")
out.append("                return Some((kba_ijk, kba_kji));")
out.append("            }")
out.append("        }")
out.append("    }")
out.append("    // Fallback: default stretch-bend params by periodic table row")
out.append("    let row_i = get_periodic_table_row(atomic_num_i);")
out.append("    let row_j = get_periodic_table_row(atomic_num_j);")
out.append("    let row_k = get_periodic_table_row(atomic_num_k);")
out.append("    for &(ri, rj, rk, kba_ijk, kba_kji) in MMFF_DFSB_TABLE.iter() {")
out.append("        if ri == row_i && rj == row_j && rk == row_k {")
out.append("            if swap {")
out.append("                return Some((kba_kji, kba_ijk));")
out.append("            } else {")
out.append("                return Some((kba_ijk, kba_kji));")
out.append("            }")
out.append("        }")
out.append("    }")
out.append("    None")
out.append("}")

with open("/Users/rliu/Projects/WebMM/src/mmff/mmff_tables.rs", "w") as f:
    f.write("\n".join(out))

print(f"Generated mmff_tables.rs with:")
print(f"  EQ levels: {len(def_data)} entries")
print(f"  Angle table: {len(angle_data)} entries")
print(f"  SB table: {len(stbn_data)} entries")
print(f"  DFSB table: {len(dfsb_data)} entries")
