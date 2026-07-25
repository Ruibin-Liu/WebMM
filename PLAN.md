# Close final 8 typing gaps — validation at 0.00% mismatch

**Status:** ✅ COMPLETE — 0% type mismatch on both validation sets (109-mol + 43-file)

## Fixes (src/mmff/mod.rs)
1. **has_c_c_neighbor / has_c_n_neighbor gated on `!n_owns_cn`**: the enamine
   (C=C) and amidine (C=N) → N_PL3 rules were catching the =N imine itself
   (its *other* neighbor has the C=C/C=N). Now excluded via `n_owns_cn`
   (this N has its own double bond to C) — the =N falls through to sp2 N_2.
   Fixed guanine atom6/8, purine atom1/7 (4 cases).
2. **N_4 charge rule gated**: quaternary-ammonium N_4 no longer fires for
   acyl/amidine/enamine Ns (charged guanidinium N → N_PL3 via amidine).
   Fixed guanidinium atom3.
3. **H_NAM gated on non-aromatic N**: aromatic-ring N-H (pyrrole/indole NPYL)
   → H_N3; H_NAM only for non-aromatic amide/amidine/enamine/sulfonamide/aniline.
   Fixed indole atom12. Broadened acyl detection to include C=C (enamine),
   fixed guanine exocyclic NH2 H → H_NAM.
4. **5-ring dicoordinate N → N5A** when adjacent to NPYL/OFUR/S_AR (else N5B).
   Fixed pyrazole atom4.

## Result
- Validation atom types: **8 → 0 (0.00%)** on 109 mols (1324 atoms); 43-file
  set also 0%. 165 tests pass, 0 warnings.
- Charges: mean |Δ| **0.0092 → 0.0003** (only 1 atom >0.15, indole borderline).
- Energy: Pearson **r 0.967 → 0.9754**, RMSD **8.65 → 7.53**.
- Caffeine minimized −123.02 (RDKit −123.49, Δ≈0.5).

Atom typing now matches RDKit 2025.09.3 exactly on the full validation set.
Remaining energy gaps are parameter-table issues (5-ring/strained-ring angles,
P/S params), not typing — see docs/validation-energy-analysis.md.
