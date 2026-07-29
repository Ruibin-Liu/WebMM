# Plan: Phase 0 — multi-seed ETKDG harness + RDKit self-consistency ceiling

## Final Status: DONE

## Context
The single-seed harness (seed 42) is too noisy (±0.005 r per molecule flip) and
r=1.0 is the wrong target — ETKDG is stochastic. Phase 0 builds the stable
measurement gate that all later phases (§5 of ETKDG_MMFF_REVIEW.md) validate
against. Measurement infrastructure only — no algorithmic change.

## Completed
- [x] **`scripts/gen_etkdg_ref.py`** — multi-seed RDKit ref (env `WEBMM_SEEDS`,
      default `42,43,100,7`); emits `{mean,stddev,seeds,n_embedded,n_atoms,embed_ok}`.
- [x] **`examples/dump_etkdg_geom.rs`** — multi-seed WebMM dump (env `WEBMM_SEEDS`,
      default `42,43,100,7`); same schema; hand-built JSON (no new serde dep).
- [x] **`scripts/validate_etkdg.py`** — compares `mean` (fallback to legacy
      `energy` for old single-seed files); reports r/RMSD/n/seed-count/outliers.
- [x] **`scripts/rdkit_self_consistency.py`** (new) — pairwise per-seed r/RMSD +
      multi-seed-mean-vs-single = the recorded ceiling.
- [x] **End-to-end run** (4 seeds, 130 mols):
  - Ceiling (RDKit vs RDKit): single-seed pairwise **r=0.990–0.996 (RMSD 4.3–6.4)**;
    multi-seed-mean vs single **r=0.996–0.998 (RMSD 2.9–4.2)** → target ≈0.997.
  - WebMM multi-seed mean: **r=0.8609, RMSD=24.66** (cf. single-seed r=0.8568 —
    the new, more stable baseline). Real gap to ceiling ≈0.13.
  - Harness now surfaces per-molecule seed-variance (e.g. DMF stddev 7.3:
    seeds give 14.5/−5.1/−0.4/1.3 → flag WebMM embedding instability for that mol).
- [x] `cargo build --release --example dump_etkdg_geom`; `cargo test` 191 pass;
      `cargo clippy --all-targets` 0 warnings. No `src/` changes → WASM unaffected.

## Notes
- `src/` untouched this turn (example + Python only); prior turn's ETKDG fixes
  (B3 shipped, D1 shipped, B1 reverted) remain in place.
- The multi-seed mean is the gate for Phases 1–6; the ceiling (~0.997) is the
  target, not 1.0.
