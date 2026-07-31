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
| Types 52, 79 | 5 | 100% | 5/5 | `scripts/val_set_new5/` |
| RDKit bulk.sdf | 32 | 32/32 | 30/32 | `scripts/val_set_bulk/` |
| **TOTAL** | **228** | **228/228** | **226/228 (99.1%)** | |

## 2 Remaining Energy Gaps (both ≤ 0.021 kcal/mol)

| Molecule | dE | Root Cause |
|----------|-----|------------|
| 3204 (tosylate-amidinium salt) | -0.021 | Salt has large canceling electrostatic terms; 3-dp charge rounding amplifies |
| 2941 (furanium cation) | +0.018 | Extremely strained O+= ring bond (dr=-0.196 Å); 3-dp param rounding amplifies |

Both are within the RDKit verbose (3-dp displayed params) vs CalcEnergy
(full-precision internal params) precision gap. For salts and extremely
strained bonds, the 3-dp rounding amplifies to ~0.02 kcal/mol. These
cannot be fixed with parameter changes — they require exact full-precision
param storage matching RDKit's internal values.

## Deep-Dive Fix History (this session)

| Compound | Original dE | Final dE | Root Cause |
|----------|------------|----------|------------|
| 1847 (azo sulfonate) | +0.67 | 0.0 | N=N and N_2-C_AR bond params used estimation (wrong r0 cascaded into SB) |
| 3204 (salt) | -0.16 | -0.021 | C_2-S_3, N_PL3-H bonds used estimation |
| 2941 (furanium) | -27.3 | +0.018 | Missing C_AR-O_2P bond, angles, AND ring torsion (37,37,51,37) |
| oxenium_formyl_fluoride | +0.55 | 0.0 | C-F bond AND F-C-H angle (5,3,11) |
| 2289 (tryptophan) | +4.9 | 0.0 | N_AM-O_3 bond used estimation |
| 2805 (phloroglucinol) | -0.44 | 0.0 | C_AR-C_AR Single r0 wrong |
| 3283 (ammonium) | +0.053 | 0.0 | C_AR-Cl bond kb wrong |
| 1424 (chloropurine) | +0.030 | 0.0 | C_AR-Cl bond kb wrong |
| 1760 (hydrazide) | +0.036 | 0.0 | N_3-N_AM, N_3-H bonds used estimation |

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).
