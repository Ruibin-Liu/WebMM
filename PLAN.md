# Improve macrocycle 4D distance-geometry embedding (root cause #3 / Phase 5)

**Status:** PLANNED (investigation-first; substantial core-DG work)

## Context — why this is the remaining blocker
The macrocycle-amide twisting chain is now fixed downstream of the embedding:
- torsion prefs are built (#1: `num_ring_bonds` counts only small rings);
- the pre-refinement hard gates are gone (Phase 2: GATE2/GATE3 advisory);
- `minimize_etkdg` now runs, with the V2=8.0 amide prefs wired in.
Yet amides still twist because **the 4D embedding hands refinement a catastrophically
bad starting geometry** (E₀ ≈ 1575–5419 kcal/mol for the 66-atom #9143 macrocycle),
and gradient descent can't escape that basin (Phase 3 proved 300→3000 iters + 3× step
still lands in E≈1000–2000 local minima). This is a **distance-geometry quality**
problem, independent of torsions/gates/refinement.

## Current embedding algorithm (`generate_initial_coords_from_bounds`)
1. **Random per-pair distance sampling**: `dist[i][j] = lo + (hi-lo)*rng()` for every
   pair, independent of all other pairs.
2. **Classical MDS**: build double-centered Gram matrix `B = -0.5·J·D²·J`,
   Jacobi-eigendecompose (100 sweeps), take top-4 eigenvectors × √eigenvalue as 4D coords.
3. **Negative-eigenvalue fallback**: dimensions with λ≤0 get `1.0 - 2.0*rng()` (uniform
   noise in [-1,1]) — injects random noise rather than a structured coordinate.
4. One 4D minimization (`minimize_4d_first`, 400 iters) → 4D collapse → project → 3D refine.

Triangle smoothing on the *bounds* IS applied (`bounds.smooth_triangle_inequality`,
line 4112), so bounds are metric-consistent. The problem is that **independent per-pair
sampling from consistent bounds still yields a non-metric distance matrix**, and at 66
atoms (~2000 pairs) the inconsistency is severe → many negative eigenvalues → lots of
noise injection → a high-energy starting conformation.

## Phase 1 — Localize the failure (instrument, then revert)
Before changing the algorithm, confirm WHERE the embedding goes wrong for the macrocycle:
- **Bounds sanity**: after `build_distance_bounds` + smoothing, are the macrocycle
  ring-closure / 1-4 / long-range bounds reasonable (not over-tight / contradictory)?
  Does `find_rings` returning the macrocycle 4× (18,18,19,19) create redundant/
  conflicting bounds? Log bound ranges for key atom pairs.
- **Sampling quality**: how many negative eigenvalues does the Gram matrix have for a
  typical macrocycle sample (vs a small molecule)? How often does the noise fallback fire?
- **Embedding vs bounds**: after the raw metric embedding (pre-minimization), what's the
  energy and how many bounds are violated? Is the 4D minimization reducing it or stuck?
This tells us whether to attack bounds, sampling, projection, or minimization.

## Phase 2 — Fix `find_rings` over-enumeration (foundational, likely prerequisite)
`find_rings` (BFS non-tree-bond ring finder, `graph.rs:333`) is NOT a true SSSR — it
returns the #9143 macrocycle **4×** (18,18,19,19) and over-counts small rings (4× size-5
where RDKit sees 5,5,6,6). This inflates ring-dependent distance bounds (ring closure,
1-3, 1-4) with redundant/conflicting entries. Replace with a proper SSSR (or at minimum
dedup identical rings + select a minimal basis). High leverage (fixes bounds for ALL
ringed molecules) but **high regression risk** — ring detection underpins aromaticity,
atom typing, torsion classification. Must validate the full 130-mol benchmark + 190 tests
after. Defer to after Phase 1 confirms bounds are actually the issue.

## Phase 3 — Improve distance sampling for metric consistency (highest leverage)
The core algorithmic weakness. Candidate approaches, ranked:
- **(a) Sample at the upper bound** (or a fixed fraction of [lo,hi]) instead of uniform
  random — produces a more consistent, "expanded" distance matrix that embeds better.
  Cheap to try; RDKit-style.
- **(b) Metropolis / coordinated sampling**: sample distances in an order that respects
  already-sampled ones (e.g., shortest-path-consistent). More work, better quality.
- **(c) Use the metric-embedding output itself** as the sample (i.e., embed from the
  smoothed *upper-bound* matrix directly, skip random sampling) — deterministic, always
  metric-consistent. This is a known DG variant worth benchmarking.
Pick based on Phase 1 evidence + macrocycle amide result.

## Phase 4 — Improve negative-eigenvalue handling
Replace the `1.0 - 2.0*rng()` noise fallback (line ~1605) with something structured:
zero-fill, or a small jitter, or project only the positive-eigenvalue subspace and place
remaining dims at 0. Noise injection is almost certainly hurting large-molecule embeddings.

## Phase 5 — Strengthen 4D minimization
If the embedding is still rough after Phases 2–4: more `minimize_4d_first` iterations
(currently 400) for macrocycles, or replace the ad-hoc gradient descent with the existing
L-BFGS optimizer (already used for MMFF) for faster convergence from bad starts. Phase 1
data (is the 4D minimizer stuck?) decides if this is needed.

## Phase 6 — Verify (hard gates)
- `examples/macrocycle_amide.rs`: all 4 amides planar (O=C–N–X ∈ ~0–20°/160–180°) across
  seeds 0xf00d/1/7/42/100/12345; `[C@@H]` stereocenter correct sign; starting E₀ dropped
  from ~5000 to a reasonable range.
- **Regression (must not change):** 130-mol benchmark — 0% type mismatch, r≈1.0000 vs
  RDKit 2026.03.4; `cargo test` 190 pass; `cargo build` + `cargo clippy --all-targets` 0 warnings.
- Optionally commit the #9143 SDF + add a `#[test]` asserting amide planarity.

## Phase 7 — Update CODE_STATUS.md

## Risks / constraints
- **Highest risk = Phase 2 (find_rings/SSSR)** — touches a foundation used by typing,
  aromaticity, torsions. Gate it on the full benchmark; only do it if Phase 1 shows bounds
  are the issue.
- Phases 3–5 change the embedding for **all** molecules, not just macrocycles — every
  change must pass the 130-mol regression (types + energy) and the 190 ETKDG tests.
- No MMFF changes; no API changes; no new dependencies.
- This is genuinely hard (RDKit's macrocycle DG is the product of years of tuning);
  set expectations: incremental improvement, possibly not full RDKit parity in one pass.
- Stop-and-assess after each phase against the benchmark — do not batch risky changes.

## Prerequisite (already in place)
#1 (`num_ring_bonds`) + Phase 2-of-prior-plan (advisory gates) are merged, so once the
embedding improves, amide planarization + correct stereocenter should follow without
further torsion/gate work.
