# Plan: Review-fix pass 2 — remove production debug prints, tighten aniline bounds guard

## Context
Review of HEAD (b3e02bd) found two issues in the review-fix pass:
1. `minimize_etkdg` (src/etkdg/mod.rs) — the core 3D ETKDG minimizer on the
   public/WASM-exposed default embedding path — emits production debug output:
   an unconditional `eprintln!("MIN n=…")` on every embedding and a `DBGLS stall`
   print on line-search stalls. Clippy does not flag `eprintln!`, so these slipped
   through; they are leftover diagnostics from earlier passes (a7a24b0, fb30419).
2. `test_aniline_hh_bounds` only asserts `H-H lower > 1.6`, which the *old buggy*
   bounds also satisfy (126° split → lower ≈1.795); the test would not catch a
   regression to the 114°/123° split. The H-H 1-3 bound center should encode ≈120°
   (measured center 1.7831 Å = exactly 120°; buggy 126° → 1.8348 Å).

## Phase 1 — Remove production debug prints (src/etkdg/mod.rs)
- Remove the `eprintln!("DBGLS stall …")` block in the line-search-failure branch
  of `minimize_etkdg`; keep the `stall`/history-reset/break logic (functional).
- Remove the unconditional `eprintln!("MIN n=…")` at the end of `minimize_etkdg`.
- Remove the now-unread `iters_done` / `converged` bindings (init, increment, and
  `converged = true;` — keep the `break`) so `cargo clippy --all-targets` stays
  at 0 warnings.
- `n` stays (used by the final coord copy-back loop). No behavior change.

## Phase 2 — Tighten `test_aniline_hh_bounds` (src/etkdg/mod.rs)
- Replace the weak `lower > 1.6` assertion with an angle check: from the printed
  H-H and N-H 1-3 bound centers, compute the implied H-N-H angle and assert
  `|angle − 120°| < 2°` (sp2 N with no double bond → flat 2π/3). The buggy
  126° split fails by 6°, so this is a real regression guard.

## Phase 3 — Docs per workflow
- Overwrite `PLAN.md` (this file).
- Update `CODE_STATUS.md`: new Recently Completed entry + Current Focus sentence.

## Verification
- `cargo test` — 199/199 (aniline tests ×3 runs, deterministic).
- `cargo clippy --all-targets` — 0 warnings.
- `python3 scripts/benchmark_mmff.py --no-speed` — 230/230 PASS (no FF change).
- `cargo test test_aniline_hh_bounds -- --nocapture` — prints bounds; assertion
  passes on the new center (120.0°) and would fail on the old split (126.0°).
- Commit the pass.

## Out of scope
- Test-module `eprintln!`/`println!` diagnostics (intentional, cfg(test)-only).
- `src/lib.rs` pre-existing TODO comments (test code, not this pass).
- ETKDG embedding quality / harness numbers (unchanged by this pass).
