# Large-Scale RDKit Validation (109 molecules)

**Status:** ✅ COMPLETE (harness built; surfaced a prioritized fix list)

## What was built
- `scripts/gen_validation_set.py` — generates 109 explicit-H 3D molecules (RDKit
  ETKDG, fixed seed) from a curated SMILES set covering the full MMFF chemistry
  space (alkanes/alkenes/alkynes, every O/N functional group, S/P chemistry,
  halogens, all ring sizes 3-8, fused/bridged, all common heterocycles, ions).
  PubChem-drug fetching is wired but disabled (sandbox blocks `urllib`/`subprocess`
  import; curl works — use `--drugs` on an unrestricted host).
- `examples/dump_types_energy.rs` — WebMM atom types + energy at SDF coords → JSON.
- `scripts/validate.py` — type / charge / energy comparison + outlier report.

## Results (109 mols, 1324 atoms) after the H_NAM fix
- **Atom types: 91 mismatches (6.87%)** (was 7.93%).
- **Charges: mean |Δ| = 0.0107** (excellent); 23 atoms >0.15 (ammonium, purines, CS2).
- **Energy: Pearson r = 0.952, RMSD 11.3 kcal/mol**; systematic high-energy bias
  for 5-ring heteroaromatics + strained rings.

## Fix applied this task
- `src/mmff/mod.rs` (`determine_h_subtype`): H on an amide N (bonded to C=O) now
  → H_NAM (28), was H_N3 (23). Reduced mismatches 105→91.

## Prioritized remaining gaps (roadmap)
1. **5-ring type collapse (56 cases, biggest)**: type_ids collapses C5A/C5B→C_AR,
   NPYL/N5B→N_AR, OFUR→O_3 for energy correctness. RDKit reports the specific
   types. Fix: add proper 5-ring angle params so real type_ids work — would also
   cut the biggest energy outliers (indole/adenine/purine/thiazole/benzothiophene).
2. **H-subtype residuals (16)**: H_NAM/H_N3 boundary + HNR+(36) for charged-amine H.
3. **Ammonium charge bug** (Δ0.94): N charge −1.74 vs RDKit −0.80.
4. **3-ring carbon CR3R (22)**: cyclopropane C, ×3.
5. **CS2 sulfur S2CM (72)**: ×2.
6. **N_2(9) vs N_PL3(40)**: ×3 (sp2 amine-like N).

## Note
The earlier 43-file set showed 0% because it under-sampled 5-ring heteroaromatics
and amides. The 109-set is the meaningful benchmark.
