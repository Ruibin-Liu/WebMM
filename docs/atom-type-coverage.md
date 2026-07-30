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
| Types 52, 79 | 5 | 100% | 3/5 | `scripts/val_set_new5/` |
| RDKit bulk.sdf | 32 | 31/32 | 24/32 | `scripts/val_set_bulk/` |
| **TOTAL** | **228** | **227/228** | **218/228** | |

## All 86 Real Types Now Covered

The final 2 types were found in this session:

| Type | Element | Description | SMILES | Molecule |
|------|---------|-------------|--------|----------|
| 52 | H | HO=+ (H on oxenium O+=) | `[OH+]=C` | Oxenium methaniminium |
| 79 | N | N5 (general N in 5-ring w/ alpha+beta N heteroatoms) | `C[n+]1cn[nH]c1` | Methyl triazolium |

Type 79 requires a charged 5-ring (e.g. triazolium) where an N at degree 2
has pyrrole-like N neighbors at both alpha and beta positions — achievable
only with a pyrrole-like N (2 pi electrons) + a pyridinium-like N+ (1 pi
electron) giving exactly 6 aromatic pi electrons.

## Remaining Energy Gaps (10 molecules)

| Molecule | dE | Root Cause |
|----------|-----|------------|
| 2289 (tryptophan deriv.) | +4.9 | Likely stretch-bend or charge param |
| 2941 (furanium cation) | -27.3 | O+= (type 51) params in ring context |
| 2805 (phloroglucinol) | -0.44 | Moderate param gap |
| 1847 (azobenzene) | -0.43 | Type mismatch: azo N=9 vs 40 |
| 3204 (sulfonamide salt) | -0.16 | Small param gap |
| methyl_triazolium | -0.76 | Charge/param fine-tuning |
| oxenium_formyl_fluoride | +0.55 | C=O+ bond param |
| 1424, 1760, 3283 | <0.06 | Rounding-level |

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).
