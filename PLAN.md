# Plan: Phase 1 — port RDKit RNG (investigated → PARKED)

## Final Status: PARKED — premise was wrong; low value for the r≈1.0 goal

## What was investigated
- **The old review doc claimed RDKit uses `std::mt19937`.** Verified against RDKit
  2025.09.3 source `Code/RDGeneral/utils.h:36`:
  `typedef boost::minstd_rand rng_type;` + `boost::uniform_real<double>`
  (`variate_generator`). From the boost 1.84 header, `boost::minstd_rand` =
  `linear_congruential_engine<uint32_t, 48271, 0, 2147483647>` — the **48271
  Park–Miller LCG** (`x = 48271·x mod 2^31−1`), mapped to `[0,1)` by boost's
  `uniform_real`. **Not MT19937.** So the roadmap's "Phase 1 = port MT19937" was
  based on the doc's false claim.
- **Bit-exact validation is not cleanly available:** `pickRandomDistMat` is not
  exposed in the RDKit Python API (only `DoTriangleSmoothing`, `EmbedBoundsMatrix`),
  and there are no local boost headers to compile a C++ reference probe.

## Why parked
- **It does not move the r metric.** RNG selects *which* conformer a given seed
  yields, not embedding *quality*. The r≈1.0 gate is the **multi-seed mean** (ceiling
  ~0.997), which is RNG-independent. My own roadmap said Phase 1 "alone won't move r."
- **Risk without validation:** a source-only port of boost's LCG + `uniform_real`
  float-generation could be subtly wrong (seeding edge cases, multi-draw float
  mapping) with no way to confirm bit-exactness here.

## Outcome
- No code changes (correctness preserved; no unvalidated RNG shipped).
- `ETKDG_MMFF_REVIEW.md` corrected: E11 + §4 table + §5 Phase 1 now state the real
  target (`boost::minstd_rand`, not MT19937) and the low-value/hard-to-validate
  status.

## Recommendation (pending user decision)
Skip Phase 1; proceed to **Phase 2 (faithful bound matrix bundle)** — that is where
real r gains are, and it is fully validatable via the multi-seed harness with no RNG
dependency. (If bit-exact RNG is still wanted later: `brew install boost` + a ~15-line
C++ probe to capture the `minstd_rand`+`uniform_real` double sequence, then port +
unit-test against that vector.)

## Non-goals
Any algorithmic change this turn.
