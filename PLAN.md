# Fix #1: 5-ring type_ids collapse (and H_NAM subtype)

**Status:** ✅ COMPLETE (2026-07-25)

## #1 — 5-ring type_ids un-collapse
`type_ids` previously collapsed 5-ring heteroaromatic types to their aromatic
base (C5A/C5B→C_AR, NPYL/N5A/N5B→N_AR, OFUR→O_3) "for energy correctness".
This made WebMM report C_AR/N_AR where RDKit reports C5A/NPYL/etc. — the single
biggest type-mismatch category (56 of 91 on the 109-mol validation set).

Investigation: the collapse was a workaround for a *different* bug (the NPYL
charge −1.03 / mis-typed 5-ring N). Once that was fixed, the collapse became
**energy-neutral** — caffeine's full energy breakdown is byte-identical with or
without it (TOT −88.39; r=0.954, RMSD 11.04 unchanged across 109 mols).

Fix: `src/mmff/mod.rs` `type_ids` now uses `mmff_type_id(at)` for every type
(no base_type collapse). Validation mismatches **91 → 51 (3.85%)**. 165 tests
pass, 0 warnings. (Earlier "un-collapse breaks caffeine" was the NPYL charge
bug, not the collapse itself.)

## Remaining (cosmetic, no energy impact)
- **C5A vs C5B alpha/beta swap (13 cases)**: RDKit keys on position relative to
  pyrrole- vs pyridine-type N, not hetero-neighbor count. Both EQ-fall to C_2,
  so zero energy impact. Documented as a known limitation, not pursued.
- H-subtype residuals (16), HNR+ charged-amine H (5), cyclopropane CR3R (3),
  water OH2 vs O_3 (3), N_2 vs N_PL3 (3), CS2 S2CM (2), H_COOH vs H_ONC (2).

## H_NAM subtype (also done this session, prior commit's validation)
H on amide N → H_NAM(28), was H_N3(23).
