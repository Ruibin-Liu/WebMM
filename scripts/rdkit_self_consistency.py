#!/usr/bin/env python3
"""RDKit-vs-RDKit seed scatter — the irreducible r ceiling on the ETKDG harness.

Reads the multi-seed rdkit_ref.json (from gen_etkdg_ref.py) and reports, over the
molecules that embedded at every seed:
  - pairwise r / RMSD between each pair of seeds (single-conformer scatter), and
  - mean(all seeds) vs each single seed (the multi-seed ceiling).

This is the *target* WebMM is measured against: ETKDG is stochastic, so even
RDKit-vs-RDKit single-seed r is ~0.99, not 1.0. validate_etkdg.py should compare
multi-seed means, which lifts the ceiling to ~0.997 and cuts per-molecule flip
noise.
"""
import json, math, sys, itertools

ref = json.load(open('scripts/etkdg_val/rdkit_ref.json'))
# seeds that appear anywhere
seed_set = sorted({int(s) for v in ref.values()
                   if isinstance(v, dict) for s in v.get('seeds', {})})
if not seed_set:
    sys.exit('no per-seed energies in rdkit_ref.json (regenerate with gen_etkdg_ref.py)')

# molecules embedded at ALL seeds (clean per-seed comparison)
common = [n for n, v in ref.items()
          if isinstance(v, dict) and v.get('embed_ok')
          and all(str(s) in v.get('seeds', {}) for s in seed_set)]
print(f"Seeds: {seed_set}")
print(f"Molecules embedded at every seed: {len(common)}")


def stats(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n)
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n)
    r = cov / (sx * sy) if sx * sy > 0 else 0.0
    rmsd = math.sqrt(sum((x - y) ** 2 for x, y in zip(xs, ys)) / n)
    return r, rmsd


def col(seed):
    return [ref[n]['seeds'][str(seed)] for n in common]


print("\nPairwise single-seed scatter (RDKit vs RDKit):")
print(f"  {'pair':<16}{'r':>9}{'RMSD':>9}")
for a, b in itertools.combinations(seed_set, 2):
    r, rmsd = stats(col(a), col(b))
    print(f"  seed {a:<4} vs {b:<4}  {r:>9.4f}{rmsd:>9.2f}")

# mean over all seeds vs each single seed
means = [sum(ref[n]['seeds'][str(s)] for s in seed_set) / len(seed_set)
         for n in common]
print("\nMulti-seed ceiling (mean over all seeds vs a single seed):")
print(f"  {'comparison':<24}{'r':>9}{'RMSD':>9}")
for s in seed_set:
    r, rmsd = stats(means, col(s))
    print(f"  mean(seeds) vs seed {s:<5}  {r:>9.4f}{rmsd:>9.2f}")

# headline: mean-vs-mean would be trivially 1.0, so the actionable ceiling for a
# single-seed *reference* is mean(seeds) vs that one seed.
print("\n=> Target for validate_etkdg.py (multi-seed mean vs multi-seed mean): "
      "ceiling ~0.997.")
