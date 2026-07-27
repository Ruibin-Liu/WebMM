# ETKDG validation harness (mirror of the MMFF94s val set)

**Status:** IN PROGRESS

## Why
The 130-mol MMFF benchmark uses RDKit-supplied geometry, so it cannot catch
WebMM **ETKDG embedding** regressions. Every candidate macrocycle-embedding fix
therefore can't be validated. Build a parallel harness: embed each molecule with
WebMM ETKDG and measure embedding **quality** vs RDKit's ETKDG, so changes to the
embedding pipeline have a regression signal.

## Design (mirrors MMFF val set exactly)
- **MMFF val set**: `gen_validation_set.py` (RDKit → SDFs + `rdkit_ref.json`) +
  `dump_types_energy.rs` (WebMM → `webmm_ref.json`) + `validate.py` (compare).
- **ETKDG val set** (this): same 3-part shape, but the WebMM side **embeds** rather
  than reading fixed geometry, and the metric is **embedding quality**, not type/energy
  at a fixed geometry.

### Metric
Both RDKit and WebMM embed each molecule independently (fixed seed) and report the
**single-point MMFF energy of their own embedded conformer** (a strain measure), plus
`embed_ok` (finite, non-degenerate coords) and `n_atoms`. A good WebMM embedding has
energy comparable to RDKit's; a bad one (e.g. macrocycles) is much higher → flagged as
an outlier. This is directly analogous to the MMFF energy comparison.

## Phase 1 — `scripts/gen_etkdg_ref.py` (RDKit)
- Reuse the curated SMILES from `gen_validation_set.py`.
- For each: AddHs → ETKDGv3 embed (randomSeed=42) → sanitize (skip kekulize) →
  set aromaticity → MMFF94s props → record embedded-conformer MMFF energy.
- Write connectivity SDF to `scripts/etkdg_val/<name>.sdf` + `scripts/etkdg_val/rdkit_ref.json`
  `{name: {energy, n_atoms, embed_ok}}`.

## Phase 2 — `examples/dump_etkdg_geom.rs` (WebMM)
- Read each `scripts/etkdg_val/*.sdf`, `parse_sdf`, run WebMM `generate_initial_coords_with_config`
  (random_seed=42, macrocycle torsions on, defaults), build `MMFFForceField`, compute
  single-point energy of the embedded coords.
- Emit `{name: {energy, n_atoms, embed_ok}}` JSON to stdout (→ `webmm_ref.json`).

## Phase 3 — `scripts/validate_etkdg.py`
- Load both refs; for common molecules compare embedded energies.
- Report: count embedded, Pearson r of (rdkit_e, webmm_e), RMSD, and the worst
  outliers (webmm_e ≫ rdkit_e = bad WebMM embedding). This is the regression baseline.

## Phase 4 — Baseline run + record
- Run all three; capture baseline numbers (r, RMSD, worst outliers — expect
  macrocycles/complex molecules to dominate the outlier list).
- Commit-ish artifacts under `scripts/etkdg_val/`. `cargo build`/`clippy`/190 tests unchanged.

## Phase 5 — Update CODE_STATUS.md
- Document the harness + baseline; note it now enables safe ETKDG regression checks
  for future embedding work (incl. the macrocycle fix).

## Constraints
- New files only (`scripts/gen_etkdg_ref.py`, `scripts/validate_etkdg.py`,
  `examples/dump_etkdg_geom.rs`, `scripts/etkdg_val/`). No production-code changes.
- `uv`/pip N/A (use system `python3` which has RDKit 2025.09.3; ETKDG geometry
  quality is version-stable enough for a regression baseline — version pinned in output).
- No new Rust deps.
