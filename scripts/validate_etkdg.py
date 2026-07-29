#!/usr/bin/env python3
"""Compare WebMM vs RDKit ETKDG-embedded-conformer MMFF energy (multi-seed mean).

Both tools embed each molecule independently at K seeds and we compare the MEAN
single-point MMFF energy across seeds. Comparing means (rather than a single seed)
cuts the stochastic noise floor: RDKit-vs-RDKit single-seed r~0.99, multi-seed-mean
r~0.997 (see rdkit_self_consistency.py). So the ceiling is ~0.997, not 1.0.

Backward compatible: falls back to the legacy single-value `energy` field if a ref
file lacks `mean` (old single-seed format).
"""
import json, math, sys

rd = json.load(open('scripts/etkdg_val/rdkit_ref.json'))
us_path = sys.argv[1] if len(sys.argv) > 1 else 'scripts/etkdg_val/webmm_ref.json'
us = json.load(open(us_path))


def val(d):
    """Return the comparable energy: prefer mean-over-seeds, else legacy energy."""
    if not isinstance(d, dict):
        return None
    return d.get('mean', d.get('energy'))


def nseeds(d):
    return len(d.get('seeds', {})) if isinstance(d, dict) else 0


common = sorted(set(rd) & set(us))
pairs = []
for name in common:
    r, u = rd[name], us[name]
    if not (isinstance(r, dict) and isinstance(u, dict)):
        continue
    if not r.get('embed_ok') or not u.get('embed_ok'):
        continue
    rv, uv = val(r), val(u)
    if rv is None or uv is None:
        continue
    pairs.append((name, rv, uv))

rd_fail = sum(1 for n in common if isinstance(rd[n], dict) and not rd[n].get('embed_ok'))
us_fail = sum(1 for n in common if isinstance(us[n], dict) and not us[n].get('embed_ok'))
print(f"Molecules in common: {len(common)}")
print(f"Embedded OK (both): {len(pairs)}")
print(f"Embed failures: rdkit={rd_fail}, webmm={us_fail}")
# report seed coverage for the first embedded molecule (representative)
for name, r, u in pairs[:1]:
    print(f"Seeds compared (first mol {name}): rdkit={nseeds(rd[name])}, "
          f"webmm={nseeds(us[name])}")

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
    print(f"\nEmbedded-energy (mean-over-seeds) Pearson r: {corr:.4f}")
    print(f"Embedded-energy (mean-over-seeds) RMSD: {rmsd:.2f} kcal/mol")
    print(f"(ceiling ~0.997 — see scripts/rdkit_self_consistency.py)")
    deltas = sorted(
        ((name, re_, we, we - re_) for name, re_, we in pairs),
        key=lambda x: -abs(x[3]),
    )
    print(f"\nWorst |ΔE| (name, rdkit, webmm, Δ=webmm-rdkit):")
    for name, re_, we, d in deltas[:20]:
        print(f"  {name:<24} rdkit={re_:>9.2f} webmm={we:>9.2f} Δ={d:+.2f}")
