# Plan: Fix discrepancies found in new validation compounds

## Status: 30/41 fixed (<0.01 kcal/mol)

## Completed
- [x] **D. H on phosphorus** — type 71 for H bonded to P
- [x] **A. Sulfonium/hypervalent S** — type 17 (S+ 3 bonds), type 18 (4+ bonds)
- [x] **E. Imine N in cumulated N=C=O** — type 9 not 10 (n_owns_cn guard)
- [x] **F. Peroxide O-O + Si-N bonds** — Si-N 4.254/1.700, O-O 4.088/1.449
- [x] **G. Cumulated double bonds** — C_1=O_2, C_1=N_2, C_1=S_2 params + angle entries
- [x] **H. 4-membered rings** — CR4R-S, CR4R-O_R bonds + C-S-C angle
- [x] **I. Oxazole** — N5A-OFUR bond params
- [x] **SF4** — type 18 (removed Sp3D interception) + S-F bond params
- [x] **Nitroso N (type 46)** — N_NITROSO variant + detection + bond params
- [x] **Cumulated N (types 53/47)** — N_2Z/N_1M variants + detection + bond params
- [x] **Nitro charge fix** — is_nitro_n requires charge > 0.5

## Remaining (11 molecules)
- [ ] **B. Isonitrile** (types 60/61) — methyl_isocyanide (+53), benzonitrile_oxide (+56)
- [ ] **C. Guanidinium cation** (types 55/57/36) — guanidinium_salt (-70)
- [ ] **methyl_nitrite** — O type wrong (32→7), needs O type detection fix
- [ ] **phosphirane** — 3-membered P ring (+9.5)
- [ ] Small param gaps: methyl_azide (+0.4), SF4 (-0.9), cumulated systems (<0.6)
