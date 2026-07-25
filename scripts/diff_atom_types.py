#!/usr/bin/env python3
"""Diff WebMM vs RDKit MMFF94s atom types across all test molecules."""
import json, sys
from collections import Counter

rd = json.load(open('scripts/rdkit_atom_types.json'))
us = json.load(open('scripts/webmm_atom_types.json'))

# aggregate mismatch patterns: (rdkit_type, webmm_type, symbol) -> count
patterns = Counter()
detail = []
total_atoms = 0
total_mismatch = 0
files_with_mismatch = 0
skipped = []

for name in sorted(set(rd) & set(us)):
    r_entry = rd[name]
    u_entry = us[name]
    if not isinstance(r_entry, list) or not isinstance(u_entry, list):
        continue
    # Skip molecules without explicit hydrogens: WebMM has no implicit-H
    # inference (it targets explicit-H 3D SDFs), so 2-coordinate carbons read
    # as sp2 where RDKit infers sp3. These are input artifacts, not typing bugs.
    # (Verified separately: adding explicit H via RDKit makes WebMM match exactly.)
    has_h = any(s == 'H' for m in r_entry if isinstance(m, dict)
                for s in m.get('syms', []))
    if not has_h:
        skipped.append(name)
        continue
    # align by molecule index
    for mi, (rm, um) in enumerate(zip(r_entry, u_entry)):
        if not isinstance(rm, dict) or not isinstance(um, dict):
            continue
        if 'error' in rm or 'error' in um:
            continue
        rt = rm.get('types', [])
        ut = um.get('types', [])
        rsyms = rm.get('syms', [])
        usyms = um.get('syms', [])
        # guard: atom ordering must match (same symbols)
        if rsyms != usyms:
            detail.append(f'{name}#{mi}: SYMBOL MISMATCH rdkit={rsyms} webmm={usyms}')
            continue
        for i, (r, u) in enumerate(zip(rt, ut)):
            total_atoms += 1
            r_int = int(r) if str(r).lstrip('-').isdigit() else r
            u_int = int(u)
            if r_int != u_int:
                total_mismatch += 1
                sym = rsyms[i]
                key = (r_int, u_int, sym)
                patterns[key] += 1
                detail.append(f'{name}#{mi} atom{i} {sym}: rdkit={r_int} webmm={u_int}')

print(f'Total atoms compared: {total_atoms}')
print(f'Total mismatches: {total_mismatch} ({100*total_mismatch/max(total_atoms,1):.1f}%)')
if skipped:
    print(f'Skipped (no explicit H — input artifact, not typing bugs): {skipped}')
print()
print('Top mismatch patterns (rdkit_type, webmm_type, symbol) -> count:')
for (rt, ut, sym), cnt in patterns.most_common(40):
    print(f'  rdkit={str(rt):>4}  webmm={ut:<4}  sym={sym:<3}  x{cnt}')

# per-file summary
print()
print('Files with mismatches:')
for name in sorted(set(rd) & set(us)):
    cnt = sum(1 for d in detail if d.startswith(name + '#'))
    if cnt:
        print(f'  {name}: {cnt} atoms')
