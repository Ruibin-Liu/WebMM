# Comprehensive Review: WebMM ETKDGv3 & MMFF94s vs. RDKit

> **Re-audited against RDKit 2025.09.3 source** (`Code/DistGeom/TriangleSmooth.cpp`,
> `DistGeomUtils.cpp`, `Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`,
> `Embedder.cpp`). The previous version of this document marked several ETKDG
> issues "✅ FIXED" that are **not** fixed in the code; this version reflects the
> actual state. MMFF94s items (M1–M8) are validated independently by the
> 130-molecule MMFF energy harness (r=1.0000, RMSD=0.170) and remain solid.

## Executive summary

- **MMFF94s (M1–M8): genuinely resolved** — validated by the fixed-geometry MMFF
  harness (r=1.0000, 0 outliers >1.0). Not re-litigated here.
- **ETKDGv3: substantially NOT RDKit-equivalent.** Three items the old doc called
  "fixed" are not (E5, E11; E12 was fixed in this audit). More importantly, the old
  doc omitted a large class of **algorithmic deviations** (D1–D12 below) in the
  minimizers, 3D force field, torsion table, and macrocycle handling. These are why
  the ETKDG embedding harness sits at **r=0.857 / RMSD=24.8** vs the
  **RDKit-vs-RDKit ceiling of r≈0.99 / RMSD≈5** (see §3).

---

## 1. ETKDGv3 issues — actual status

### E5: Triangle smoothing — **NOT FIXED**
`smooth_triangle_inequality` (`src/etkdg/mod.rs:206`) has **two** deviations from
RDKit `triangleSmoothBounds`:
1. **Wrong lower-bound formula:** uses `|lower[i][k] − lower[k][j]|` instead of
   RDKit's `max(L_ik − U_kj, L_kj − U_ik)` (`diffLikUjk`/`diffLjkUik`). The current
   formula can produce infeasible (over-tight) lower bounds.
2. **Fixed-point vs single sweep:** loops `while changed` to convergence, whereas
   RDKit does a single Floyd–Warshall k-sweep (upper bound = shortest path converges
   in one sweep; lower bound is deliberately left partial).

Fixing (1) alone, or (1)+(2), **regressed** the harness (r→0.852 / 0.850) because
WebMM's downstream pipeline is co-adapted to the old bounds. Reverted; needs the
bundled migration in §5 Phase 2.

### E11: RNG — **NOT FIXED**
`Rng` (`src/etkdg/mod.rs:59`) is **Xoshiro256\*\*** seeded by splitmix64
(`rotl(s[1]*5,7)*9` output, `s[1]<<17`/`rotl(s[3],45)` update, splitmix64 constants).
RDKit uses `std::mt19937`. Bit-exact conformer reproducibility is impossible until
this is ported. (Only affects reproducibility, not average quality.)

### E12: `is_larger_sp2_atom` ring check — **FIXED (this audit)**
`src/etkdg/mod.rs:623` now requires `numAtomRings > 0`, matching RDKit's
`isLargerSP2Atom`. (Previously doubled `DIST13_TOL` for all sp² atoms with Z>13,
including non-ring ones.)

### E9: Force constants — **premise was WRONG; constants actually match**
The old note claimed "RDKit uses weight=1.0". RDKit uses
`KNOWN_DIST_FORCE_CONSTANT=100` for 1-2/1-3 and `boundsMatForceScaling*10 = 10` for
long-range (`DistGeomUtils.cpp`). WebMM's `K_12=100` / `LONG_RANGE_FORCE=10` **match
these exactly.** The real 3D-force-field deviation is D4 (see §2), not the constants.

### Other E-items (claimed fixed, status as of this audit)
- **E1/E2 (UFF electronegativity/radii):** formula structure verified; values
  match the standard UFF set. ✅
- **E3/E4 (chiral R/S, planarity check):** present in code; not re-derived in this
  audit. *Needs confirmation.*
- **E6 (ET version):** consistent — RDKit `ETKDGv3` uses `ETversion=2` (the "v3"
  comes from `useMacrocycleTorsions`/`useMacrocycle14config`); WebMM default
  `et_version=2`. ✅
- **E7/E8 (torsion ordering, bridged-ring exclusion):** *not re-verified* in this
  audit — the torsion system is a hand-coded approximation (D10), so "ordering" is
  moot until the real SMARTS table is ported.
- **E10 (timeout units):** uses seconds. ✅

---

## 2. Algorithmic deviations NOT in the old doc (D1–D12)

These are the real reasons ETKDG embedding quality trails RDKit.

| ID | Area | WebMM | RDKit |
|----|------|-------|-------|
| D1 | `set_ring_angle` SP3D2 | fixed this audit (90°) | `_setRingAngle` 90° |
| D2 | `check_and_set` | tightens (intersect) | `_checkAndSetBounds` widens (union) |
| D3 | sp² 1-3 angle | 123°/114° heuristic | flat 120° / `(2π−taken)/remaining` |
| D4 | 3D long-range term | `BASIN_THRESH`+`MAX_UPPER` filter; double-counts 1-2/1-3 | constrains **all** non-1-2/1-3/1-4 pairs, no filter |
| D5 | 1-2/1-3 pin target | bounds-matrix `[lo,hi]` | current distance ± `KNOWN_DIST_TOL` |
| D6 | planarity | harmonic impropers + extra ring/exocyclic torsion terms | UFF `InversionContrib` Fourier only |
| D7 | 4D minimizer | fixed-step descent, single pass | BFGS in `while(needMore)` |
| D8 | `computeInitialCoords` | never fails | fails on bad eigenvalues → retry |
| D9 | tetra/chiral pre-checks | advisory (no retry) | hard retry gates |
| D10 | torsion preferences | ~12 hand-coded cases | full SMARTS tables (~hundreds) |
| D11 | post-processing | torsion-snap, bond-snap, H-trilateration, H-relax, aniline-NH₂ | none (not needed) |
| D12 | misc | "conjugated"≈non-Single; no in-ring 1-4 stereo; approximate macrocycle amide +0.1 | computed `getIsConjugated()`; stereo Z/E; dedicated macrocycle SMARTS |

**What already matches RDKit well:** all 17 distance constants; weight staging
(first min `1.0/0.1`, fourth-dim `0.2/1.0`); iteration counts `400/200/300`;
attempt count `10·nAtoms`; `optimizerForceTol 1e-3`; random-distance sampling
`lb+r·(ub−lb)`; torsion-pref Fourier form; metric-embedding double-centering
(different form, mathematically equivalent); UFF bond-length formula; `basinThresh=5`.

---

## 3. The r≈1.0 target — and why it's ~0.99, not 1.0

ETKDG is stochastic. Measured **RDKit-vs-RDKit** seed scatter on the 130-mol set
(single conformer per molecule, ETKDGv3, MMFF94s single-point energy):

| pair | r | RMSD |
|------|---|------|
| seed 42 vs 43 | 0.9908 | 6.41 |
| seed 42 vs 100 | 0.9927 | 5.66 |
| seed 42 vs 7 | 0.9952 | 4.59 |
| mean(42,43,100,7) vs 42 | **0.9973** | **3.42** |

Implications:
- The **irreducible single-seed ceiling is r≈0.99 / RMSD≈5.** Chasing 1.0000 on a
  single-seed harness is chasing noise.
- WebMM at **r=0.857 / RMSD=24.8** has a **real ~0.13 algorithmic gap** to close
  (≈5× the strain it should have). This is achievable via the migration in §5.
- A **multi-seed harness** (mean over ≥4 seeds) is the correct gate: it raises the
  ceiling to ~0.997 and cuts per-molecule flip noise (currently ±0.005).

---

## 4. Corrected fix-status table

| Issue | Status |
|-------|--------|
| E1 Electronegativity | ✅ Verified |
| E2 Hybridization radii | ✅ Verified |
| E3 Chiral R/S | ⚠️ Present, not re-verified |
| E4 Planarity check | ⚠️ Present, not re-verified |
| E5 Triangle smoothing | ❌ **Not fixed** (formula + single-sweep; reverted, needs bundle) |
| E6 ET version | ✅ Consistent (ETKDGv3 = ETversion 2) |
| E7 Torsion ordering | ⚠️ Moot until D10 (hand-coded patterns) |
| E8 Bridged ring exclusion | ⚠️ Not re-verified |
| E9 Force constants | ✅ Constants match (old premise wrong; real issue is D4) |
| E10 Timeout units | ✅ Seconds |
| E11 RNG | ❌ **Not fixed** (still Xoshiro256**, not MT19937) |
| E12 is_larger_sp2_atom | ✅ **Fixed this audit** |
| M1–M8 MMFF94s | ✅ Validated by MMFF harness (r=1.0000) |

---

## 5. Road to r≈1.0 — systematic plan

**Core principle:** WebMM's pipeline is *co-adapted* to its own deviations — that's
why one-off RDKit-faithfulness fixes regress (B1 lesson). The path is a coordinated
**inside-out migration to RDKit-faithfulness**, bundled per phase, gated by a
**multi-seed harness**, with the ad-hoc workarounds (D11) feature-flagged for
removal once parity is shown.

**Phase 0 — Measurement (do first).** Convert the harness to mean-over-K-seeds
(≥4). Track the RDKit-vs-RDKit ceiling (≈0.99) as the target, not 1.0. Cheap,
removes the noise that's been masking signal.

**Phase 1 — RNG foundation (E11/B2).** Port `std::mt19937` + match `pickRandomDistMat`
draw order (`i=1..n, j=0..i`). Necessary for reproducibility; alone it won't move r
much but it unblocks Phases 2–4 validation.

**Phase 2 — Faithful bound matrix.** Bundle (land together, single commit):
B1 smoothing formula + single sweep; D2 `check_and_set` widen; D3 sp² 120°/
distribute; D12 conjugation/stereo/in-ring details. Re-tune nothing else yet —
expect a dip, recovered in Phases 3–4.

**Phase 3 — Faithful minimizers.** D7: generalize the existing L-BFGS to the 4D
kernels, run `while(needMore)` like RDKit. D8: restore `computeInitialCoords`
failure/retry. D9: restore tetrahedral/chiral retry gates. **Then** delete the D11
workarounds (torsion-snap, bond-snap, H-trilateration, H-relax, aniline) — they
exist only to compensate for the weak minimizer and should become redundant.

**Phase 4 — Faithful 3D force field.** D4: long-range over all non-(1-2/1-3/1-4)
pairs, no `BASIN_THRESH`/`MAX_UPPER` filter, no 1-2/1-3 double-count. D5: pin 1-2/1-3
to current distance ± tol. D6: replace harmonic+extra planarity terms with UFF
`InversionContrib` Fourier (forceScaling 10).

**Phase 5 — Full torsion preference table (D10).** Port `torsionPreferences.in` +
`torsionPreferences_macrocycles.in` SMARTS tables. **Largest single quality lever**
for floppy/conjugated systems. Requires a SMARTS matcher in WebMM.

**Phase 6 — Macrocycles.** Proper SSSR `find_rings`; macrocycle 1-4 bounds; the
#9143 ring-amide fix. Closes the known macrocycle gap (currently worse than RDKit
pre-fix).

**Sequencing rule:** Phases 2+3+4 are one logical migration — do not merge them
independently or intermediate states will regress (the B1 result). Land them behind
a cfg flag, validate the whole bundle on the multi-seed harness, then flip.
