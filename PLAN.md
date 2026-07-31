# Plan: ETKDG default-path outliers — small-ring cluster (start: cyclopropane)

## Goal
Lift default-path r from 0.8603 toward the ceiling by fixing the concrete outliers.
Worst cluster (stable MMFF baseline): strained 3–7-membered rings (cyclopropane +64,
aziridine +59, cyclopentene +59, cyclohexene +50, cycloheptane +87, THF +47) and
hypervalent-S compounds (sulfone/sulfonate/sulfonamide/sulfuric_acid +46–93).

## Approach (default path — the shipped pipeline; faithful path parked behind flag)
Debug one outlier deeply (geometry vs RDKit + per-term MMFF breakdown), fix the root
cause, measure on the deterministic harness, then extend to the cluster.

## Steps
- [ ] **Diagnose cyclopropane:** embed with WebMM (seed 42, default) vs RDKit; compare
      coords (bonds/angles/H) AND per-term MMFF energy (bond/angle/torsion/vdw/elec/
      sb/oop) to localize the strain (is it the C–C bonds, C–C–C angles, H–C–H, H
      out-of-plane, or 1-4/VDW?).
- [ ] **Fix root cause** in `src/etkdg/*` (likely a bounds/torsion/H-placement detail).
- [ ] **Measure:** `validate_etkdg.py` — cyclopropane Δ→0 and r up; 191 tests; clippy.
- [ ] **Extend** to aziridine / cyclopentene / cyclohexene / cycloheptane / THF.
- [ ] Then the hypervalent-S cluster.

## Non-goals
Faithful path (parked behind `rdkit_faithful`). MMFF (fixed, 228/228). No refactors.

## Gate
Default path byte-stable (0.8603 → up); 191 tests; clippy 0 (mine); deterministic.
