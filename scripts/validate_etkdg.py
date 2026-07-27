#!/usr/bin/env python3
"""Compare WebMM vs RDKit ETKDG-embedded-conformer MMFF energy — regression baseline.

Both tools embed each molecule independently (seed 42) and we compare the single-point
MMFF energy of their embedded conformers. A good WebMM embedding has energy comparable
to RDKit's; a much higher WebMM energy flags a bad embedding (strained/clashed/non-planar).
This is the ETKDG analog of scripts/validate.py (which compares at fixed RDKit geometry).
"""
import json, math, sys

rd = json.load(open('scripts/etkdg_val/rdkit_ref.json'))
us_path = sys.argv[1] if len(sys.argv) > 1 else 'scripts/etkdg_val/webmm_ref.json'
us = json.load(open(us_path))

common = sorted(set(rd) & set(us))
pairs = []
for name in common:
    r, u = rd[name], us[name]
    if not isinstance(r, dict) or not isinstance(u, dict):
        continue
    if not r.get('embed_ok') or not u.get('embed_ok'):
        continue
    if r.get('energy') is None or u.get('energy') is None:
        continue
    pairs.append((name, r['energy'], u['energy']))

rd_fail = sum(1 for n in common if isinstance(rd[n], dict) and not rd[n].get('embed_ok'))
us_fail = sum(1 for n in common if isinstance(us[n], dict) and not us[n].get('embed_ok'))
print(f"Molecules in common: {len(common)}")
print(f"Embedded OK (both): {len(pairs)}")
print(f"Embed failures: rdkit={rd_fail}, webmm={us_fail}")

if pairs:
    xs = [p[1] for p in pairs]
    ys = [p[2] for p in pairs]
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n)
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n)
    corr = cov / (sx * sy) if sx * sy > 0 else 0.0
    rmsd = math.sqrt(sum((x - y) ** 2 for x, y in zip(xs, ys)) / n)
    print(f"\nEmbedded-energy Pearson r: {corr:.4f}")
    print(f"Embedded-energy RMSD: {rmsd:.2f} kcal/mol")
    # Worst outliers — WebMM energy far above RDKit's = bad embedding.
    deltas = sorted(
        ((name, re_, we, we - re_) for name, re_, we in pairs),
        key=lambda x: -abs(x[3]),
    )
    print(f"\nWorst |ΔE| (name, rdkit, webmm, Δ=webmm-rdkit):")
    for name, re_, we, d in deltas[:20]:
        print(f"  {name:<24} rdkit={re_:>9.2f} webmm={we:>9.2f} Δ={d:+.2f}")
