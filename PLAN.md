# Fix aziridine 3-ring torsion gap

**Status:** ✅ COMPLETE

## Root cause
`find_torsions` created degenerate torsions where atom1 == atom4 in 3-membered
rings (the ring wraps around to the same atom: C2-C0-C1-C2). RDKit filters
these; WebMM didn't, inflating torsion energy.

## Fix
`src/molecule/graph.rs`: added `l != k` filter in find_torsions.
`src/lib.rs`: added `no_degenerate_torsions_in_3_ring` regression test.

## Result
- cyclopropane: Δ +0.71 → 0.00 (exact)
- aziridine: Δ +1.24 → +0.41
- 130-set: r=1.0000, RMSD 0.281 → 0.248, 0 outliers >1.0
- 190 tests pass, 0 warnings; WASM rebuilt
