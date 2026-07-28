# Plan: Fix discrepancies found in new validation compounds

## Final Status: 40/41 types match, 31/41 energy <0.01

## Completed (all atom typing patterns fixed)
- [x] **A. Sulfonium/hypervalent S** — type 17 (S+), type 18 (4+ bonds)
- [x] **B. Cumulated N** — types 53 (=N=), 47 (NIM-)
- [x] **C. Guanidinium cation** — types 57 (CGD+), 55 (NCN+), 36 (HNC+)
- [x] **D. H on phosphorus** — type 71
- [x] **E. Imine N** — type 9 (n_owns_cn guard)
- [x] **Nitroso N** — type 46
- [x] **Isonitrile** — types 60 (CID), 61 (NID)
- [x] **Nitro charge guard** — is_nitro_n requires charge > 0.5

## Completed (param fixes)
- [x] Si-N, O-O, CR4R-S, CR4R-O_R, N5A-OFUR, S-F bonds
- [x] C_1=O_2/N_2/S_2 cumulated double bonds
- [x] Cumulated angle entries (1,9,4), (4,9,27), (9,4,9), (9,4,16)
- [x] C-S-C 4-ring angle, N_2-N_2Z bond, CID≡NID bond

## Remaining (10 molecules, all param-level)
- [ ] guanidinium_salt (-74): charge computation for cation types 55/57
- [ ] benzonitrile_oxide (-29): needs type 35 for nitrile oxide O⁻ (1 type bug)
- [ ] phosphirane (+9.5): 3-membered P ring params
- [ ] methyl_nitrite (+4.3): nitroso bond/angle params
- [ ] SF4 (-0.9), cumulated systems (<0.6): small param gaps
