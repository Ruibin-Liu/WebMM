# Plan: Review-fix pass (flaky ETKDG test, GBSA robustness, repo hygiene)

## Context
Review of HEAD found: (1) `etkdg::tests::test_aniline_geometry` is flaky — aniline NH2
H-N-H embeds at 122–126° (RDKit ~117.5°) and the 110–125° assertion fails ~30% of seeds;
root cause is the sp2 1-3 angle heuristic giving the H-H pair the leftover 126° when the
center has three single bonds (aniline N in Kekulé form), plus the default random seed.
(2) GBSA has no d≈0 guard in the Born-radius passes, the born_inv floor clamp is not
reflected in `chain`, the NVE/FD test guards are much looser than the measured values,
and the LCPO containment-boundary force kink is undocumented. (3) PLAN.md is stale,
`Cargo.lock` version bump is uncommitted, and AGENTS.md describes a FastAPI/React stack
that does not match this Rust/WASM repository.

## Phase 1 — Fix flaky aniline ETKDG geometry test
- `src/etkdg/mod.rs` `build_distance_bounds`: in the sp2 `visited_centers >= 1` branch,
  when the center has NO double bond (all-single neighbors, e.g. aniline N), distribute
  the remaining angle equally (`remaining / remaining_pairs`) instead of applying the
  carbonyl 114°/123° split — matches RDKit set13Bounds (flat 2π/3 for non-ring sp2).
  Carbonyl-like centers (has_double) keep the existing split.
- Pin `random_seed` in `test_aniline_geometry` (and `test_aniline_hh_bounds`) so the
  embedding is deterministic.
- Verify: `examples/diag_aniline_seeds.rs` (temporary) — all 40 seeds pass; full
  `cargo test` twice in a row (no flakes); `cargo clippy --all-targets` 0 warnings.

## Phase 2 — GBSA robustness (`src/solvation/mod.rs`)
- Guard `d ≈ 0` in the ψ (pass 1) and born-radius-gradient (pass 3) pair loops so
  coincident atoms cannot produce NaN.
- Make `chain` consistent with the `born_inv` floor clamp: when the clamp is active the
  derivative is 0 (piecewise chain).
- Tighten `gb_nve_stays_finite` assertion (10% → measured 0.07% with margin) and the FD
  gradient threshold to the measured max error with margin.
- Document the LCPO S_ij containment-boundary force kink in a code comment and in
  CODE_STATUS.md known limitations (model-fidelity, not smoothed).

## Phase 3 — Repo hygiene
- Overwrite `PLAN.md` (this file) per workflow.
- Update `CODE_STATUS.md` (Recently Completed + Current Focus + Known Risks/Issues).
- Align `AGENTS.md` stack overview with the actual Rust/WASM toolchain (cargo,
  wasm-pack; python3 + RDKit for validation scripts); keep the mandatory
  CODE_STATUS.md/PLAN.md workflow and scope constraints.
- Commit: Cargo.lock 0.5.0→0.6.0 bump + all fix changes.

## Verification
- `cargo test` ×2 — all pass, no flake.
- `cargo clippy --all-targets` — 0 warnings.
- `python3 scripts/benchmark_mmff.py --no-speed` — 230/230 PASS (no MMFF change).
- `examples/diag_aniline_seeds.rs` — 40/40 seeds within bounds (temporary, removed).

## Out of scope
- GBSA WASM API (native-only by design, not part of this pass).
- Smoothing the LCPO model (deviates from Weiser-Tsui-Case; documented instead).
- Other ETKDG embedding quality (4D embedding, macrocycles) — separate efforts.
