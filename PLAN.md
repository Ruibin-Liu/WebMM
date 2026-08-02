# Plan: Decouple MMFF atom typing from aromaticity + take over ETKDG aromaticity fix

## Context
The other session's `graph.rs` WIP broadens the N aromatic-candidate rule from 5-rings
to any ring (`ring_bonds == 2` without the `ring.len() == 5` gate). This is intended to
improve ETKDG planarity for xanthine/purines, but it regresses 15 molecules because MMFF
atom typing uses the aromaticity flag to override carbonyl-C and amide-N-H typing:
- Carbonyl C: C_2(3) → C_AR(37)  →  neighbor N_AM charge −0.49 → −0.671 (35 kcal divergence)
- Amide N-H:  H_NAM(28) → H_N3(23)

**Root cause:** MMFF atom typing is coupled to the aromaticity flag. RDKit's MMFF typer
types ring carbonyl carbons as C_2 and amide N-H as H_NAM **regardless of aromaticity**.
At HEAD this is latent (graph.rs doesn't over-mark), but it breaks under the new
aromaticity perception.

## Phase 1 — Harden MMFF atom typing (decouple from aromaticity flag)

### Fix A: Carbonyl C → C_2, not C_AR
`src/mmff/mod.rs` carbon cascade, before `(6, _, true, _) => C_AR`:
```rust
// Carbonyl C (C=O double bond) → C_2 even in aromatic rings.
// RDKit types ring carbonyl C (xanthine, uracil, 2-quinolone) as C_2, not C_AR.
(6, _, true, _) if double_bond_partners.contains(&8) => MMFFAtomType::C_2,
```

### Fix B: Amide/amidine N-H → H_NAM regardless of aromaticity
`src/mmff/mod.rs` `determine_h_subtype`: split the `!n_is_aromatic && (acyl || ...)` gate
so the acyl/amidine case (N with a C neighbor that is C=O/C=N/C=C) yields H_NAM even when
the N is aromatic. Pyrrole (aromatic, no acyl) still → H_N3; sulfonamide/aniline unchanged.

### Verification
- `cargo test` — all pass (no MMFF regressions; fix is no-op at HEAD since graph.rs
  doesn't over-mark carbonyl C as aromatic there)
- `python3 scripts/benchmark_mmff.py` — 230/230 still match RDKit <0.01 kcal/mol
- With graph.rs WIP applied: xanthine types revert to HEAD values (C_2, H_NAM, charge −0.49)

## Phase 2 — Evaluate graph.rs aromaticity change for ETKDG

After Phase 1 makes MMFF typing robust, re-apply the graph.rs N-aromaticity broadening
and measure ETKDG harness impact:
- Does xanthine's ETKDG quality improve (was +14.5)?
- Do the 15 regressed molecules recover (MMFF energy now correct)?
- If net-positive and no ETKDG regressions → keep; else → refine the aromaticity rule

## Out of scope
- Metadynamics WASM API + site (already done this session)
- Further ETKDG phases (4D embedding, torsion SMARTS) — separate efforts
