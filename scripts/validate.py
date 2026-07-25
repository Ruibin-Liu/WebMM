#!/usr/bin/env python3
"""Compare WebMM vs RDKit on the validation set: atom types, energies, charges."""
import json, sys, math
from collections import Counter

rd = json.load(open('scripts/val_set/rdkit_ref.json'))
us = json.load(open('scripts/val_set/webmm_ref.json'))

# ---- atom types ----
total_atoms = 0
type_mismatch = 0
type_patterns = Counter()
type_detail = []
charge_l1 = 0.0
charge_n = 0
charge_big = []

# ---- energy ----
energy_pairs = []  # (name, rdkit, webmm)

common = sorted(set(rd) & set(us))
for name in common:
    r = rd[name]
    u = us[name]
    if not isinstance(r, dict) or not isinstance(u, dict):
        continue
    if 'error' in r or 'error' in u:
        continue
    rt, ut = r['types'], u['types']
    rc, uc = r.get('charges', []), u.get('charges', [])
    if len(rt) != len(ut):
        continue
    mm = 0
    for i, (a, b) in enumerate(zip(rt, ut)):
        total_atoms += 1
        if int(a) != int(b):
            type_mismatch += 1
            mm += 1
            type_patterns[(a, b)] += 1
            if len(type_detail) < 30:
                type_detail.append(f"{name} atom{i}: rdkit={a} webmm={b}")
    # charges
    for cr, cw in zip(rc, uc):
        d = abs(cr - cw)
        charge_l1 += d
        charge_n += 1
        if d > 0.15:
            charge_big.append((name, cr, cw, d))
    energy_pairs.append((name, r['energy'], u['energy']))

print(f"Molecules compared: {len(energy_pairs)}")
print(f"Atoms compared: {total_atoms}")
print(f"Type mismatches: {type_mismatch} ({100*type_mismatch/max(total_atoms,1):.2f}%)")
if type_patterns:
    print("  Top type-mismatch patterns (rdkit, webmm) -> count:")
    for (a, b), c in type_patterns.most_common(15):
        print(f"    rdkit={a} webmm={b}: {c}")
    print("  detail:")
    for d in type_detail:
        print(f"    {d}")

# charges
if charge_n:
    print(f"\nCharge mean |Δ|: {charge_l1/charge_n:.4f} over {charge_n} atoms")
    charge_big.sort(key=lambda x: -x[3])
    print(f"  atoms with |Δcharge| > 0.15: {len(charge_big)} (showing top 12)")
    for name, cr, cw, d in charge_big[:12]:
        print(f"    {name}: rdkit={cr} webmm={cw} Δ={d:.3f}")

# energy
if energy_pairs:
    print(f"\nEnergy comparison ({len(energy_pairs)} mols):")
    # correlation
    xs = [p[1] for p in energy_pairs]
    ys = [p[2] for p in energy_pairs]
    n = len(xs)
    mx = sum(xs)/n; my = sum(ys)/n
    cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))/n
    sx = math.sqrt(sum((x-mx)**2 for x in xs)/n)
    sy = math.sqrt(sum((y-my)**2 for y in ys)/n)
    corr = cov/(sx*sy) if sx*sy > 0 else 0
    rmsd = math.sqrt(sum((x-y)**2 for x,y in zip(xs,ys))/n)
    print(f"  Pearson r: {corr:.4f}")
    print(f"  RMSD: {rmsd:.2f} kcal/mol")
    # worst outliers (by abs delta)
    deltas = sorted(((name, r, w, w-r) for name, r, w in energy_pairs), key=lambda x: -abs(x[3]))
    print(f"  Worst energy deltas (name, rdkit, webmm, Δ):")
    for name, r, w, d in deltas[:15]:
        print(f"    {name:<24} rdkit={r:>9.2f} webmm={w:>9.2f} Δ={d:+.2f}")
