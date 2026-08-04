# Plan: Live animation pacing — make the motion actually visible

## Context
User: "I don't see the compound move at all" on the live MD/metad animation.
Measured via CDP: the viewer DID move (atom 0 displaced ~3.5 Å over the run,
0.63 Å max per frame) but the whole animation completed in ~1.15 s headless
(~0.45 s at 60 fps rAF) — over before the eye registered it.

## Phase 1 — Pacing fix (site/index.html, runMD + runMetad)
- Frame count: `chunk = max(10, ceil(nSteps / 120))` → ~120 frames per run
  (was ~25), ≥10 fs per frame so per-frame motion stays visible.
- Pace: `setTimeout(1000 / fps)` with `fps` from the existing speed selector
  (default 30) instead of `requestAnimationFrame` (which runs at display
  refresh and flashes the run past).
- Expected: default 2000-step MD → ~119 frames ≈ 4 s with ~0.1–0.2 Å/frame.

## Phase 2 — Docs per workflow
- Overwrite `PLAN.md`; update `CODE_STATUS.md`.

## Verification (CDP headless)
- MD live: 119 frames, per-frame max displacement ~0.17 Å, final
  E=−101.2/T=384K (matches batch), playback enabled, zero console errors.
- MetaD live: 119 frames, 40 hills, FES span 4.4 (matches batch).
- Note: headless SwiftShader renders ~40–70 ms/frame → ~8.7 s wall in
  headless; on a GPU browser the pacing yields ~4 s.
- Commit, push.

## Out of scope
- Physics (dt, step count) — unchanged; determinism preserved.
- Any engine changes.
