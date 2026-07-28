# Add conjugated-diene (C=C–C=C) torsion preference

**Status:** IN PROGRESS

## Bug
`match_torsion_pattern` has a comment "C=C-C=C (conjugated diene) → trans" but its
condition is `is_double_bond(a2,a3) && sp2(h2) && sp2(h3)` — that matches the
**cumulated** case (central bond double, C=C=C), not a conjugated diene (central
bond **single**). So butadiene/isoprene's C=C–C=C central bond gets **no torsion
preference** and embeds ~90° twisted (low MMFF energy only because MMFF barely
penalizes the twist — butadiene E≈12 at 87° vs RDKit ≈6 planar). RDKit embeds it
s-trans (~180°).

## Phase 1 — Add the correct conjugated-diene pattern (DONE)
In `match_torsion_pattern`: for a central bond a2–a3 that is **single**, with both
a2,a3 = sp² C and non-aromatic (i.e. a conjugated C=C–C=C / C=C–C=O / styrene-side
single bond), return a trans-favoring + planarizing Fourier term:
`([1,-1,1,1,1,1], [2.0, 4.0, 0,0,0,0])` → k=1 favors trans (180°), k=2 enforces
planar (0/180). Min at trans; cis at +4; 90° at +10 kcal/mol. Also fix the
mislabeled comment on the existing cumulated pattern.

## Phase 2 — Verify
- `examples/diag_butadiene.rs`: dihedral C0-C1-C2-C3 → ~180° (s-trans) across
  seeds; energy → ~6 (RDKit-like).
- ETKDG harness: r ≥ 0.61, RMSD ≤ 47 (no regression); diene outliers (isoprene)
  should drop.
- `cargo test`: 191 pass, 0 ignored; `cargo build` + `cargo clippy`: 0 warnings.

## Phase 3 — Update CODE_STATUS.md

## Constraints
- One new pattern + a comment fix in `match_torsion_pattern`. No other changes.
- If the sp²-sp²-single condition over-constrains some conjugated system (styrene/
  enone), narrow the condition rather than weakening globally — watch the harness.
