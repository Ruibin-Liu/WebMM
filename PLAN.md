# Plan: Phase 2+3+4 coordinated bundle — RDKit-faithful ETKDG migration

## Status
- **M0 (scaffolding) — DONE.** `rdkit_faithful: bool` added to `ETKDGConfig`;
  `generate_initial_coords_with_config` now dispatches to `_default` (unchanged
  body) vs `_rdkit` (passthrough). Dump example toggles via `WEBMM_RDKIT_FAITHFUL`.
  Verified: default r=0.8603, rdkit r=0.8603, **0 molecule diffs** (true passthrough);
  191 tests; my code clippy-clean (1 pre-existing warning is in `src/mmff`, other
  session). WASM rebuilt. Next: **M1** (empirical A/B — start with D4 long-range).

## Goal
Lift the ETKDG embedding harness from **r=0.8603 / RMSD=24.75** toward the
**~0.997 ceiling** (RDKit-vs-RDKit multi-seed mean; single-seed ceiling ~0.99).
Real gap to close: ~0.13. Target: r ≥ 0.95, RMSD ≤ ~8.

## Why a coordinated bundle (not incremental)
This session proved **Phase 2 (bounds) and Phase 3 (4D minimizer) each REGRESS in
isolation** (Phase 2: catastrophic on P=O / sulfones; Phase 3: r→0.85 across all
ablations). Root cause: the 3D stage (`minimize_etkdg`) + the ad-hoc workarounds
(torsion-snap, bond-snap, H-trilateration, H-relax) are **co-adapted** to the
current under-converged bounds+minimizer. There is no single-component win.
⇒ The bounds, minimizer, and 3D force field must migrate **together**.

## Architecture: runtime flag + parallel functions (shipped r never at risk)
- Add `pub rdkit_faithful: bool` (+ per-concern bools for granular A/B) to
  `ETKDGConfig` (default `false`).
- `generate_initial_coords_with_config` dispatches: if `rdkit_faithful`, call a new
  `generate_initial_coords_rdkit(mol, config)`; else the **unchanged** current path.
- New RDKit-faithful functions live **alongside** the current ones (suffix `_r`),
  not interleaved — so the working pipeline is untouched until the flip.
- The dump example reads `WEBMM_RDKIT_FAITHFUL` (and granular flags) from env, so
  old-vs-new is measurable in **one build** on the deterministic harness.
- **Only when the new path beats the old on the harness does it become default and
  the old code get deleted.**

## Milestones (each a measurable checkpoint; gate = deterministic multi-seed r vs 0.8603)

**M0 — Scaffolding (no behavior change).** Add the flag(s) + dispatch + env toggle +
a trivial `_r` passthrough that just calls the current functions. Verify r unchanged.
~1–2 h.

**M1 — Empirical A/B search (lead with the un-tested, cheapest levers).** Implement
each concern behind its flag and measure it ALONE (old bounds/4D/3D elsewhere):
1. **D4 alone** (3D long-range: drop `BASIN_THRESH`/`MAX_UPPER` filter; stop
   double-counting 1-2/1-3) — *cheapest, possibly high-impact (currently
   under-constrains loose pairs).*
2. **D5 alone** (pin 1-2/1-3 to current distance ± tol, not bounds `[lo,hi]`).
3. **D6 alone** (UFF `InversionContrib` Fourier planarity replacing harmonic
   impropers + extra ring/exocyclic torsion terms) — *biggest single piece.*
4. **Phase 4 = D4+D5+D6** combined.
5. Re-confirm **Phase 2** (B1+D2+D3+**D14**+D12-conj) and **Phase 3** (D7) alone on
   the deterministic base (prior numbers were partly on a noisy/older base).
~6–10 h (D6 dominates).

**M2 — Combination search.** From M1's per-concern deltas, try the promising
combos: **D4+D5** (cheap 3D-FF subset), **Phase 3+4**, **Phase 2+3+4 (full)**. Find
the **minimal set that lifts r above 0.8603**. If a subset wins, prefer it.
~2–4 h.

**M3 — Polish the winner.** Whatever combo wins: fix the **latent D14 bug**
(`visited_centers` per-pair → snapshot `in_ring` once/atom; required for D3) ,
restore **D8** (`computeInitialCoords` failure/retry) + **D9** (tetra/chiral retry
gates) in the new path, then strip now-redundant **D11** workarounds from the new
path. Re-measure after each. ~3–5 h.

**M4 — Ship.** If new-path r ≥ target: flip `rdkit_faithful` default to `true`,
delete the old parallel functions + workarounds, WASM rebuild, final harness.
~1–2 h. If new-path r < target but > old: optional partial ship behind a kept flag.

## Component inventory (effort / risk)
| Item | Concern | Effort | Risk |
|------|---------|--------|------|
| D4 | 3D long-range: all non-(1-2/1-3/1-4) pairs, no filter, no double-count | low | low |
| D5 | 1-2/1-3 pin to current dist ± tol | low | low |
| D6 | UFF `InversionContrib` (C0+C1·sinW+C2·cos2W, forceScaling 10) | **high** | **high** (per-atom params, grad) |
| B1 | smoothing formula `max(L_ik−U_kj,…)` + single Floyd–Warshall sweep | low | low |
| D2 | `check_and_set` widen (≡ `set_lower`/`set_upper`) | low | low |
| D3 | sp² angles (ring distribute / non-ring flat 120) | low | needs D14 |
| D14 | `in_ring` snapshot once/atom (latent bug) | low | low |
| D12-conj | `getIsConjugated()` for 5-ring squish | medium | approximable |
| D7 | 4D L-BFGS via shared `lbfgs_minimize` | low (helper designed) | low |
| D8/D9 | retry gates | medium | medium (changes accept rate) |
| D11 | remove torsion-snap/bond-snap/H-trilat/H-relax/aniline | low | medium (tuning) |

## Risks & mitigations
- **Co-adaptation persists even for the bundle** → mitigate via M1/M2 search; if the
  full bundle still underperforms, fall back to **targeted outlier fixes** on the
  current pipeline (the prior "grind," now with a deterministic harness).
- **D6 (UFF inversion) is the long pole** → implement D4/D5 first (cheap, may
  already help); treat D6 as its own sub-task with a finite-diff gradient check.
- **Determinism** → the aromatic-sort fix is uncommitted; **commit it first** so the
  trustworthy harness isn't lost mid-migration.
- **Macrocycle regressions** → keep the macrocycle-specific code on the old path
  initially (gate the new path on `!has_macrocycles`); fold in via Phase 6 later.

## Effort estimate
~15–25 h focused work. M1 (esp. D6) is the bulk. The empirical A/B search (M1/M2) is
designed to find the **smallest** winning change and may cut this substantially.

## Stop / decision points
- After M1: if **no single concern** improves r, the bundle (M2 full) is mandatory.
- After M2: if **even the full bundle** < 0.8603, abandon RDKit-faithfulness; pivot
  to targeted outlier fixes on the custom pipeline (Option B).
- Ship only when new-path r ≥ target on the deterministic harness, 191 tests green,
  0 clippy.

## Non-goals (this turn)
Implementation. This is scoping only. Phase 1 (RNG) stays parked. Macrocycle
amides (Phase 6) deferred.
