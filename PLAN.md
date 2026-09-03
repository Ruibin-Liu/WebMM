# Plan: GFN-FF review fixes round 4 — 3-ring angle branch + ring torsions

## Context (DONE)

Code review of the working tree flagged the spiro angle condition
(`ringsj+ringsk.eq.102`) as a typo — but diffing against the upstream xtb
source (grimme-lab/xtb gfnff_ini/ini2/eg, fetched fresh) proved the `102` is
genuine: xtb's `ringsatom` returns the SENTINEL 99 for atoms in no ring, so
102 = 99 + 3 means "3-ring center, one 3-ring neighbor, one acyclic neighbor".
Our `ring_size` uses 0 for acyclic, so the ported branch was dead. Same review
session exposed (via a new methylcyclopropane vs-xtb test) that the torsion
setup had only xtb's acyclic branch — the whole `lring` ring-torsion rule set
(fr3/fr4/fr5/fr6, cis minima, terminal-flank case, CB7 fix, .not.lring guard
on the extra gauche torsion, sp3-specials ordering) was missing.

## Fixes

- mod.rs: map acyclic ring_size 0 → xtb sentinel 99 at the angle-spiro check.
- mod.rs: full ring-torsion branch; helpers ringsbond() / ring_common()
  (ringstors + ringstorl semantics) over the existing rings_all data; sp3
  specials moved after both branches (xtb order); .not.lring guard added.
- detect_bonds hypercoordination filter: review FALSE ALARM — matches xtb
  icase-2 getnb literally (per-PAIR skip, no distance trimming); comment now
  documents that xtb's final topology list is the full nbf and the two
  coincide for normal organic valences.
- Test fixture: pyr_h2o.xyz moved from /tmp into src/gfnff/ (include_str!).
- Cleanups: dead rabd0 decl, duplicated nh/no lets, needless hb_ab clone,
  ipis Floyd-Warshall hoisted out of the per-pi-system loop, stale eg2_rnr
  comment.

## Verification (DONE)

- New tests_ring3::methylcyclopropane_vs_xtb: all seven terms exact
  (torsion +0.007939443 vs xtb +0.007939436401; was +0.0335 — 4x too large;
  angle +0.045553381 vs +0.045553364714) + FD 3.7e-6 kcal/mol/Å.
- cargo test 231/231, clippy 0, wasm32-unknown-unknown check clean.

## Out of scope

- Metals/eta, PBC (excluded by request); non-zero-eg1 xtb single point.
- Overcoordinated (>4/6 candidate) topologies: single-list approximation
  documented; xtb's transient icase-2 list not replicated.
