#!/usr/bin/env python3
"""Regenerate TOR94_TABLE / TOR94S_TABLE in torsion.rs from the complete extracted TSVs.
The Rust table stores (tor_type, i, j, k, l, V1, V2, V3) sorted by (tor_type, j, k, i, l)."""
import re

def load_tsv(path):
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 8:
            continue
        try:
            tt, it, jt, kt, lt = int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4])
            v1, v2, v3 = float(p[5]), float(p[6]), float(p[7])
        except ValueError:
            continue
        rows.append((tt, it, jt, kt, lt, v1, v2, v3))
    # sort by (tor_type, j, k, i, l) to match tor_binary_search's entry_key
    rows.sort(key=lambda r: (r[0], r[2], r[3], r[1], r[4]))
    return rows

def emit(rows):
    out = []
    for (tt, it, jt, kt, lt, v1, v2, v3) in rows:
        def f(x):
            s = repr(x)
            return s
        out.append(f"    ({tt}, {it}, {jt}, {kt}, {lt}, {f(v1)}, {f(v2)}, {f(v3)}),")
    return "\n".join(out)

tor94 = load_tsv('mmff_params_extracted/defaultMMFFTor.tsv')
tor94s = load_tsv('mmff_params_extracted/defaultMMFFsTor.tsv')
print(f"# TOR94: {len(tor94)} entries (i=5 count: {sum(1 for r in tor94 if r[1]==5)})")
print(f"# TOR94S: {len(tor94s)} entries (i=5 count: {sum(1 for r in tor94s if r[1]==5)})")

src = open('src/mmff/torsion.rs').read()
def replace_table(src, name, rows):
    pat = re.compile(r'(const ' + name + r': &\[(?:u8, )+u8, f64, f64, f64\] = &\[\n)(.*?)(\n\];)', re.DOTALL)
    new = pat.sub(lambda m: m.group(1) + emit(rows) + m.group(3), src, count=1)
    return new

src = replace_table(src, 'TOR94_TABLE', tor94)
src = replace_table(src, 'TOR94S_TABLE', tor94s)
open('src/mmff/torsion.rs', 'w').write(src)
print("wrote src/mmff/torsion.rs")
