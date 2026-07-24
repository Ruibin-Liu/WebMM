# Fix MMFF Formal Charge Distribution to Match RDKit

**Status:** ✅ COMPLETE (2026-07-23) — all 163 tests pass, 0 clippy errors, WASM build verified

**Goal:** Replace WebMM's naive uniform residual charge spreading with RDKit's two-stage
MMFF formal charge model (Halgren MMFF.V equation 15), so partial charges for carboxylates,
sulfonates, N-oxides, nitro groups, guanidinium, and phosphates match RDKit 2025.09.3.

## What was done
- Implemented `compute_mmff_formal_charges` following RDKit's per-type formal charge rules
- Replaced `calculate_bci_charges` with equation 15: `pChg = (1-M·v)·q0 + v·sumFormal + sumBci`
- Fixed `mmff_bond_type` (Aromatic → bt=0, matching RDKit `getMMFFBondType`)
- Verified partial charges match RDKit exactly for acetate, N-oxide, nitromethane, sulfonate, ammonium
