# Plan: Fix playback race — traj.nFrames on null (line 494 Uncaught TypeError)

## Context
User report: `WebMM/:494 Uncaught TypeError: Cannot read properties of null
(reading 'nFrames')` in `togglePlay`'s interval callback. Sequence: ▶ Play
starts `playTimer`; then Optimize/Embed/loadMol runs `traj = null;
disablePlayback()` — but `disablePlayback()` did not stop the interval or
reset `playing`, so the next tick read `traj.nFrames` on null. Uncaught
because interval callbacks are outside any try/catch.

## Phase 1 — Fix (site/index.html)
- `disablePlayback()`: clear the interval, reset `playing`, restore the ▶
  icon (was left as ❚❚), then disable buttons/slider as before.
- Interval body: guard `if (!traj || traj.nFrames < 1) { disablePlayback();
  return; }` (defense in depth).

## Phase 2 — Docs per workflow
- Overwrite `PLAN.md` (this file).
- Update `CODE_STATUS.md`.

## Verification (CDP headless Chrome)
- MD → ▶ → Optimize while playing: zero console errors; optimize completes.
- MD → ▶ → loadMol: zero errors; play button disabled with ▶ icon.
- Commit, push (Pages redeploys).

## Out of scope
- Any other UI polish; Rust/WASM changes (site-only fix).
