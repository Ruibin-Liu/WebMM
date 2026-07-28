# ETKDG torsion barrier-crossing (torsion-snap basin-hop)

**Status:** ✅ COMPLETE — torsion-snap basin-hop added; diene pref re-enabled; butadiene planar; harness r 0.6105 → 0.6375, RMSD 47.5 → 45.2

## Problem
L-BFGS gets stuck at torsional barriers from the 4D-projection start (butadiene
seed 42 → ~102°, can't cross to 180°; this is *why* the conjugated-diene torsion
pref was reverted — it pushed butadiene INTO the barrier and raised its energy).
max_attempts re-embeds don't help: the metric embedding of a linear diene
consistently starts ~87°, so every restart converges to the same ~100° basin.

## Fix (Phase 1)
Targeted **torsion-snap basin-hop** in `generate_initial_coords_with_config`,
after `minimize_etkdg`: for the torsion prefs most twisted from their preferred
minimum, rotate the appropriate fragment around the central bond to snap the
dihedral toward a preferred value, then re-minimize; keep only if total energy
drops. Directly crosses torsional barriers (unlike a coordinate shake, which
doesn't specifically rotate around a bond).

- New helper `rotate_fragment_around_bond(coords, mol, j, k, angle)`: BFS the
  k-side fragment (atoms reachable from k without crossing j), rotate them around
  the j→k axis by `angle` (Rodrigues).
- After minimize_etkdg: for each torsion pref (i,j,k,l), if the dihedral is
  >~45° from the nearer of {0,180} (the planar/trans minima most prefs want),
  snap-rotate the k-fragment to that minimum, re-run minimize_etkdg; accept only
  if `etkdg_energy` decreases (line-search-style guard).

## Phase 2 — Verify
- `examples/diag_butadiene.rs`: butadiene seed 42 dihedral → ~180° (s-trans),
  energy → ~6 (RDKit-like); consistent across seeds.
- ETKDG harness: r ≥ 0.61, RMSD ≤ 47 (no regression); ideally small gain.
- Re-enable the conjugated-diene torsion pref (Phase 1 of the prior plan) since
  barriers can now be crossed — confirm butadiene planar + harness ≥ 0.61.
- `cargo test`: 191 pass; `cargo build` + `cargo clippy`: 0 warnings.

## Phase 3 — Update CODE_STATUS.md

## Constraints
- Accept-only-if-better guard → can't regress a molecule (worst case: no-op).
- Bounded cost: a few snap+reminimize cycles per attempt (minimize_etkdg is the
  expensive part; keep the count small).
- If it doesn't cross the butadiene barrier or regresses the harness, revert and
  document (barrier-crossing is genuinely hard; RDKit uses years of DG tuning).
