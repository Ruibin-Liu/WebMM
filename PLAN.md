# Fix remaining typing gaps + document energy outliers

**Status:** ✅ COMPLETE — validation types 51→8 (0.60%), energy r=0.954→0.967

## Energy-outlier documentation
Added `docs/validation-energy-analysis.md`: categorizes the 15 worst energy
outliers (indole/adamantane/thiazole/adenine/etc.) by family (5-ring
heteroaromatics, strained rings, P/S chemistry, guanidinium charge) with
root-cause hypotheses and fix paths.

## Typing fixes (src/mmff/mod.rs, params.rs, atom_types.rs)
1. **5-ring C5A/C5B alpha/beta**: restructured to two-pass — classify N/O/S
   first, then C by neighbor type (C adjacent to NPYL/OFUR/S_AR or flanked by
   two heteroatoms → C5A; else C5B). Fixed 13.
2. **C5A-vs-C_AR / NPYL-vs-N_AR etc. (the un-collapse, prior commit)** — already done.
3. **CR3R (MMFF 22)**: sp3 C in 3-membered ring. New type. Fixed 3 (cyclopropane).
4. **HNRP (MMFF 36)**: H on sp3 ammonium N+. New type. Fixed 5 (ammonium,
   trimethylammonium).
5. **S2CM (MMFF 72)**: S in C(=S)=S. New type + rule. Fixed 2 (CS2).
6. **HOS (MMFF 33)**: H on O-S (sulfonic acid). New type. Fixed 1.
7. **OH2 bug**: water-O rule fired for O-S/P-H hydroxyls; now requires both
   neighbors be H. Fixed 3.
8. **H_COOH for P-OH**: H on O-P → H_COOH (RDKit reuses it). Fixed 2.
9. **H_NAM/H_N3 rule corrected** (was backwards): H_NAM for amide/amidine/
   sulfonamide/aniline N-H; H_N3 for aromatic-ring N-H (pyrrole) and simple
   amines. Added amidine (C=N) and sulfonamide (SO2) detection.
10. **Amidine/guanidinium N → N_PL3**: added has_c_n_neighbor. + HNRP ordering
    (amidine charged-N H → H_NAM, not HNR+).

## Result
- Validation atom types: **51 → 8 (0.60%)**. 165 tests pass, 0 warnings.
- Energy: Pearson **r 0.954 → 0.967**, RMSD **11.04 → 8.65**.
- Caffeine minimized **−115.5 → −123.02** (RDKit −123.49, Δ≈0.5).

## Remaining 8 (exotic RDKit MMFF nuances, minimal energy impact)
- guanidinium neutral amidine N → N_PL3 (1; WebMM gives N_4).
- guanine/purine 6-ring N → N_2(9) (4; RDKit uses N_2 not N_AR for these; needs
  WebMM aromaticity perception + RDKit's exact ring-N rule).
- guanine exocyclic NH2 H → H_NAM (2).
- pyrazole N5A vs N5B (1).
