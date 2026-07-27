# Fix remaining 7 energy outliers — target RMSD < 0.3 kcal/mol

**Status:** ✅ COMPLETE — 6/7 fixed, RMSD 0.441→0.251 kcal/mol

## Results
- 6 of 7 outliers fixed (all Δ < 0.3 except purine):
  - dimethyl_disulfide, CCl4, pyridazine, furan, pyrazole, iodobenzene → Δ ≈ 0
  - trimethylphosphine → Δ 0.0, methylphosphonic_acid → Δ +0.2
- **Remaining: purine (+1.12)** — stretch-bend sb_type classification bug
  (compute_stretch_bend_type returns sb_type 2 for mixed Double/Single bond
  angles where RDKit uses sb_type 1, causing kba to use DFSB default 0.3
  instead of specific STBN table entries 0.61/0.227). Needs RDKit source
  comparison for the exact getMMFFStretchBendType formula.

## Fixes applied
1. `src/mmff/bond.rs`: 13 RDKit-verbose-verified bond param corrections/additions
   (S-S, C-S, C-Cl, N-N aromatic, C-P single×2, P=O, C5A-OFUR, NPYL-N5A,
   N5A-C5B, C-I, C=N double, C-N single)
2. `src/mmff/mmff_tables.rs`: DFSB stretch-bend row canonicalization fix
   (was skipping asymmetric-row angles like P-O-H entirely)

## Metrics
- Pearson r: **1.0000** (rounds to 1.0)
- RMSD: **0.251 kcal/mol** (was 0.441)
- Outliers >1.0: **1** (purine, was 7)
- 108/109 molecules within 1.0 kcal/mol of RDKit
- Atom types 0.00% mismatch; charges mean|Δ|=0.0003
- 165 tests pass, 0 warnings; WASM rebuilt
