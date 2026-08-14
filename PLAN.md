# Plan: Playground v1.1 — bond grabbing & dihedral twisting

## Context
User feedback on the live page: atoms drag fine, but bonds can't be dragged or
even selected. The natural bond interaction for a physics toy is a **dihedral
twist**: grab a bond, drag perpendicular to it, and the smaller side of the
molecule rotates around the bond axis (methyl spins, conformations change, PE
chart responds). Ring bonds must be excluded — twisting them tears the ring.
Pure site-side feature, NO engine changes: rotation is applied through the
existing `MDLive::set_atom_position` export; strain safety reuses the existing
+100 kcal/mol auto-release valve.

## Phase 1 — site/playground.html: bond pick + twist
- Parse the SDF bond block on load → `bonds[]` + adjacency; mark ring bonds
  via BFS: removing bond (i,j) still leaves i connected to j ⇒ ring bond.
- Picking: pointerdown with no atom within 28 px → nearest NON-RING bond by
  2D point-to-segment distance, threshold 14 px (atoms always win).
- Twist: signed perpendicular pixel delta (relative to the bond's projected
  2D direction) → incremental Δθ = 0.01 rad/px, clamped ±15°/tick; rotate the
  smaller BFS side around the current bond axis (Rodrigues) via one
  `set_atom_position` per side atom; same path while paused (pointermove).
  Incremental rotation tracks the cursor without drift/runaway.
- Visual: while a bond is grabbed, the rotating side's atoms render in the
  accent color (per-frame colorfunc rebuild already in place).
- Strain safety: capture `peAtGrab` at grab; shared strainCheck auto-release
  applies to bond grabs too. Grabbing stays disabled during metad runs.
- `grabbed` becomes a tagged union {type:'atom'|'bond'}; pgState reports it;
  add test hooks pgBonds() / pgBondMidScreen(i,j) / pgDihedral(a,b,c,d).
- Hint/status copy updated ("drag an atom · twist a bond").

## Verification (CDP headless, extend /tmp/pg_drive.mjs)
- Press at a non-ring bond midpoint (ethanol C–C / caffeine N–CH3): grabs the
  bond; drag twists — endpoint atoms stay <0.05 Å, side atoms move >0.3 Å,
  the bond's dihedral changes, PE stays finite.
- Ring-bond midpoint press (caffeine ring): grabs nothing.
- Hard yank on a bond: strain auto-release fires.
- Orbit still works on empty space; zero console errors.
- Gates: cargo test 205/205, clippy 0 (no Rust change expected).

## Docs per workflow
- This file overwrites PLAN.md. After completion: CODE_STATUS.md entry +
  Current Focus update; README playground blurb mentions bond twisting.

## Out of scope
- Bond stretching, hover pre-highlight, beat-the-optimizer scoring.
- Spring-based dragging (engine restraint wrapper) — separate v1.5 item.
- Engine/WASM changes; MetaDLive perturbation.
