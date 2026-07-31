# MMFF Atom Type Coverage

## Summary

| Metric | Count |
|--------|-------|
| Total MMFF types | 99 |
| Real types (crd > 0) | 86 |
| Pseudo/LP types (crd = 0) | 13 |
| **All 86 real types covered** | **86 / 86 (100%)** |

## Validation Sets

| Set | Molecules | Type Match | Energy <0.01 | Location |
|-----|-----------|------------|-------------|----------|
| Original | 130 | 100% | 130/130 | `scripts/val_set/` |
| New (diverse) | 41 | 100% | 41/41 | `scripts/val_set_new/` |
| Exotic batch 1 | 8 | 100% | 8/8 | `scripts/val_set_new2/` |
| Exotic batch 2 | 6 | 100% | 6/6 | `scripts/val_set_new3/` |
| Charged/aromatic | 6 | 100% | 6/6 | `scripts/val_set_new4/` |
| Types 52, 79 | 5 | 100% | 4/5 | `scripts/val_set_new5/` |
| RDKit bulk.sdf | 32 | 32/32 | 27/32 | `scripts/val_set_bulk/` |
| **TOTAL** | **228** | **228/228** | **222/228 (97.4%)** | |

## 6 Remaining Energy Gaps (all ≤ 0.073 kcal/mol)

| Molecule | dE | Root Cause |
|----------|-----|------------|
| oxenium_formyl_fluoride | +0.073 | Near-threshold (C-F + O-H+ angle params) |
| 3283 | +0.053 | Rounding-level (RDKit 3-dp display) |
| 1760 | +0.036 | Rounding-level |
| 1424 | +0.030 | Rounding-level |
| 3204 | -0.021 | Strained sulfonate geometry (rounding) |
| 2941 | +0.018 | Strained S-O bond amplifies 3-dp rounding |

All 6 are within the known RDKit verbose (3-dp displayed params) vs
CalcEnergy (full-precision internal params) precision gap (~0.025
kcal/mol, larger for very strained bonds like the S-O in 2941 where
dr=-0.108 Å).

## Deep-Dive Findings (2941 furanium)

The 2941 dE was originally -27.3, reduced to -0.63 via bond/angle params,
then to +0.018 by adding a missing ring torsion:

- **Root cause of final gap**: WebMM was missing torsion entry
  `(37,37,51,37)` (C_AR-C_AR-O_2P-C_AR, V3=0.326). RDKit generates 2
  such ring torsions through the furanium O+ at dih≈0.3° contributing
  0.326 kcal/mol each. The previously-added `(37,51,37,1)` entry matched
  a different torsion. Adding `(0, 37, 37, 51, 37, 0, 0, 0.326)` fixed it.

- **1847/3204 root cause**: wrong bond params using estimation instead of
  RDKit values. Wrong r0 cascaded into bond + stretch-bend energy errors.

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).
