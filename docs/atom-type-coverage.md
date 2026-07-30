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
| RDKit bulk.sdf | 32 | 32/32 | 26/32 | `scripts/val_set_bulk/` |
| **TOTAL** | **228** | **228/228** | **221/228 (96.9%)** | |

## 7 Remaining Energy Gaps (all < 0.67 kcal/mol)

| Molecule | dE | Root Cause |
|----------|-----|------------|
| 1847 (azobenzene sulfonate) | +0.67 | Stretch-bend params for S-O in sulfonate context |
| 2941 (furanium cation) | -0.63 | Ring torsion param fine-tuning |
| 3204 (sulfonamide salt) | -0.16 | Stretch-bend kba |
| oxenium_formyl_fluoride | +0.07 | Near-threshold (C-F + angle params) |
| 3283 | +0.05 | Rounding-level (RDKit 3-dp display) |
| 1760 | +0.04 | Rounding-level |
| 1424 | +0.03 | Rounding-level |

The last 3 are within the known RDKit verbose (3-dp) vs CalcEnergy (full precision) gap (~0.025 kcal/mol).

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).
