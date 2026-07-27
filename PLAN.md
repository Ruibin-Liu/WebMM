# Fix remaining 7 energy outliers — target RMSD < 0.3 kcal/mol

**Status:** ✅ COMPLETE — ALL 7 fixed, RMSD 0.441→0.238 kcal/mol, 0 outliers >1.0

## Results
- ALL 7 outliers fixed; 109/109 molecules within 1.0 kcal/mol of RDKit
- **0 outliers >1.0 kcal/mol** (was 7)
- Worst remaining: ammonium (+0.99), caffeine (+0.95) — both <1.0

## Fixes applied
1. `src/mmff/bond.rs`: 13 RDKit-verbose-verified bond param corrections/additions
   (S-S, C-S, C-Cl, N-N aromatic, C-P single×2, P=O, C5A-OFUR, NPYL-N5A,
   N5A-C5B, C-I, C=N double, C-N single)
2. `src/mmff/mmff_tables.rs`: DFSB stretch-bend row canonicalization fix
   (was skipping asymmetric-row angles like P-O-H entirely)
3. `src/mmff/stretch_bend.rs`: sb_type canonicalization — swap bt_ij/bt_jk
   when type_i > type_k so bond_type_1 is the lower-type peripheral bond
   (matching RDKit's getMMFFStretchBendType canonical ordering). Fixed purine.

## Metrics
- Pearson r: **1.0000**
- RMSD: **0.238 kcal/mol** (was 0.441)
- Outliers >1.0: **0** (was 7)
- max|Δ|: **0.99** (was 3.87)
- Atom types 0.00% mismatch; charges mean|Δ|=0.0003
- 165 tests pass, 0 warnings; WASM rebuilt
