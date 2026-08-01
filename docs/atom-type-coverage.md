# MMFF Atom Type Coverage

## Summary

| Metric | Count |
|--------|-------|
| Total MMFF types | 99 |
| Real types (crd > 0) | 86 |
| Pseudo/LP types (crd = 0) | 13 |
| **All 86 real types covered** | **86 / 86 (100%)** |
| **Molecules matching RDKit <0.01 kcal/mol** | **230 / 230 (100%)** |

## Validation Sets

| Set | Molecules | Type Match | Energy <0.01 |
|-----|-----------|------------|-------------|
| Original | 130 | 100% | 130/130 |
| New (diverse) | 41 | 100% | 41/41 |
| Exotic batch 1 | 8 | 100% | 8/8 |
| Exotic batch 2 | 6 | 100% | 6/6 |
| Charged/aromatic | 6 | 100% | 6/6 |
| Types 52, 79 | 5 | 100% | 5/5 |
| RDKit bulk.sdf | 32 | 100% | 32/32 |
| Aryl phosphines (P–C_AR) | 2 | 100% | 2/2 |
| **TOTAL** | **230** | **230/230** | **230/230 (100%)** |

## Final Fix History (this round)

| Compound | dE | Root Cause |
|----------|-----|------------|
| 2941 (furanium) | +0.018 → 0 | S_3-C_AR bond used estimation (4.200/1.779) vs RDKit (3.565/1.765) |
| 3204 (salt) | -0.021 → 0 | **Inter-fragment electrostatic bug**: WebMM computed electrostatics for disconnected pairs; RDKit only does intra-fragment. Fixed via fragment BFS + excluded_pairs. |
| 3283 (ammonium) | +0.053 → 0 | C_AR-Cl bond kb wrong (3.5→3.378) |
| 1424 (chloropurine) | +0.030 → 0 | C_AR-Cl bond kb wrong |
| 1760 (hydrazide) | +0.036 → 0 | N_3-N_AM, N_3-H bonds used estimation |
| oxenium_formyl_fluoride | +0.073 → 0 | F-C-H angle (5,3,11) fallback wrong |

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).
