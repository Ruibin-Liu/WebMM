# Plan: Fix demo viewer rendering after ETKDG embed (V2000 column shift)

## Context
User report: after clicking "Embed 3D (ETKDG)" the molecule viewer shows small
dots instead of ball-stick. Root-caused via the CDP/headless-Chrome harness
(local http server over pkg/ + SwiftShader WebGL):

`renderCoords` rebuilds the SDF atom lines as
`'    ' + x(10) + y(10) + z(10) + original.substring(31)` — the 4-space prefix
plus substring(31) shifts every fixed-width V2000 field +4 columns. 3Dmol
parses by fixed columns, so it read garbage: elements "215"/"081"/"753" and
coordinates x=2, y=9660, z=5058 (from the misaligned decimal points) → atoms
render as tiny default-colored dots, mostly off-screen. This affects ALL
flows that call renderCoords (Embed, Optimize, MD), not just Embed.

## Phase 1 — Fix renderCoords (site/index.html)
- V2000 fixed columns: x[0,10) y[10,20) z[20,30) element[30,33) (0-indexed).
- Build the line as `${x}${y}${z}${lines[4+i].substring(30)}` — no leading
  spaces, splice at index 30 (element column), preserving the rest of the
  atom record (charge/mass/stereo fields).

## Phase 2 — Docs per workflow
- Overwrite `PLAN.md` (this file).
- Update `CODE_STATUS.md` (Recently Completed).

## Verification (CDP headless Chrome)
- After embed/optimize/MD, viewer model elements = {C:8, N:4, O:2, H:10}
  (caffeine) with correct coordinates; zero console errors.
- `cargo test`/clippy untouched (site-only change).
- Commit, push (Pages redeploys).

## Out of scope
- Any Rust/WASM changes (rendering was a site-side formatting bug).
- Other viewer polish.
