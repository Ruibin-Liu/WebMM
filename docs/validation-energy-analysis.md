# Validation: Energy Outliers vs RDKit

Benchmark: **109 molecules** (`scripts/val_set/`), each a RDKit ETKDG-generated
explicit-H 3D structure (seed 42). Energies compared at the **same geometry**
(WebMM `calculate_energy` vs RDKit `MMFF.CalcEnergy`), isolating force-field
parameters from the optimizer/ETKDG.

**Overall (after 2026-07-25 fixes):** Pearson r = 0.987, RMSD 5.48 kcal/mol
(was r=0.954, RMSD 11.3). 5-ring aromatics within ~1-3 kcal/mol; P/S bond
parameter gaps fixed.

### Fixed: 5-ring aromatic angle inflation (biggest win)
Root cause: `perceive_aromatic_bonds` ran once, so fused systems whose Kekulé
form gives the 5-ring only 4 pi electrons (indole) failed aromaticity on the
first pass and kept Kekulé bonds. That made the MMFF angle_type 1/2 instead of
0, selecting wrong 120-128° MMFFANG entries. Fix: iterate perception to a
fixpoint. Energy r 0.975 -> 0.984, RMSD 7.53 -> 6.05.

### Fixed: P_4 / C=S bond parameters
P_4 (phosphine oxide / phosphonate) bonds were hitting the estimator (only P_3
entries existed). Added RDKit-extracted values (C_3-P_4, P_4-O_3, P_4-O_2/
O_CO2) and corrected C_2-S_2 (thioacetone C=S). methylphosphonic_acid,
triethyl_phosphate, trimethylphosphine_oxide, thioacetone dropped out of the
worst list. r 0.984 -> 0.987, RMSD 6.05 -> 5.48.

## Remaining worst outliers
| Molecule | RDKit | WebMM | Δ | Cause |
|----------|------:|------:|---:|-------|
| adamantane | 36.49 | 75.97 | +39.5 | cage: torsion=35 + vdw=32 — RDKit per-term not exposed; same params/geometry give different totals (needs RDKit-source cross-ref) |
| norbornane | 32.32 | 46.02 | +13.7 | strained bridged (torsion+vdw) |
| caffeine | -107.27 | -95.14 | +12.1 | residual purine |
| carbon_disulfide | 17.30 | 6.32 | -11.0 | S=C=S (cumulated double bonds; finite-diff unreliable) |
| nicotine | 43.87 | 54.01 | +10.1 | residual |
| theophylline | -110.61 | -101.12 | +9.5 | residual purine |

These remaining gaps need RDKit-source-level per-term comparison (RDKit's Python
API doesn't expose the MMFF term breakdown). Adamantane's torsion=35 with the
same params/geometry/torsion-count as RDKit yet different total is the open
question.

## Worst offenders (|ΔE| > 12 kcal/mol)

| Molecule | RDKit | WebMM | Δ | Family |
|----------|------:|------:|---:|--------|
| indole | 21.61 | 62.17 | **+40.6** | 5-ring heteroaromatic (fused) |
| adamantane | 36.49 | 75.97 | **+39.5** | strained cage (C-C-C ~60°) |
| methylphosphonic_acid | −101.29 | −63.83 | +37.5 | P chemistry |
| guanidinium | −120.99 | −158.38 | **−37.4** | charged (charge bug) |
| thiazole | 57.92 | 95.23 | +37.3 | 5-ring heteroaromatic |
| adenine | 27.33 | 62.38 | +35.1 | purine (fused 5+6) |
| purine | 27.33 | 57.41 | +30.1 | purine |
| benzothiophene | 37.10 | 62.33 | +25.2 | 5-ring heteroaromatic (fused) |
| carbon_disulfide | 17.30 | 38.25 | +21.0 | S=C=S (exotic) |
| caffeine | −107.27 | −88.30 | +19.0 | purine |
| triethyl_phosphate | −51.37 | −33.92 | +17.5 | P chemistry |
| imidazole | 3.69 | 19.04 | +15.4 | 5-ring heteroaromatic |
| methanesulfonamide | −62.02 | −77.18 | −15.2 | sulfonamide |
| norbornane | 32.32 | 46.02 | +13.7 | strained bridged |
| theophylline | −110.61 | −97.70 | +12.9 | purine |

## Root-cause families

### 1. 5-membered heteroaromatics (indole, thiazole, adenine, purine, benzothiophene, imidazole, caffeine, theophylline) — WebMM systematically TOO HIGH
All Δ positive (+13 to +41). These molecules' rings are typed with the correct
specific IDs (C5A/C5B/NPYL/N5B/OFUR) after the un-collapse, but their EQ-level
fallback (C5A→C_2, N5B→N_2, etc.) lands on WebMM's **sp2 (120°) angle params**
where the rings need their true ~108° interior angles. The angle term is the
dominant contributor. **Fix path:** add explicit 5-ring angle parameters
(MMFFANG entries for C5A/C5B/NPYL/N5B/OFUR type triples) or extend the
equivalence levels so 5-ring types fall to C_AR/N_AR instead of C_2/N_2.

### 2. Strained rings (adamantane, norbornane, cyclopropane, cyclobutane) — TOO HIGH
Compressed C-C-C angles (~60–90°). WebMM applies tetrahedral 109.47° params
everywhere (MMFF CR3R/CR4R types for 3/4-ring carbons are unimplemented for
3-rings, and the ring-angle overrides are incomplete). **Fix path:** add CR3R
typing (MMFF 22) + ring-size angle overrides.

### 3. Phosphorus (methylphosphonic_acid, triethyl_phosphate) — TOO HIGH (less negative)
P=O / P-OH / P-OR params. WebMM types P as P_3(26) but the surrounding O/C bond
and angle parameters are sparse. **Fix path:** P-specific bond/angle params.

### 4. Sulfur exotic (carbon_disulfide, methanesulfonamide) — mixed
CS₂ is S=C=S (RDKit type S2CM/72, WebMM S_2/16) — wrong type + sparse params.
Methanesulfonamide is too negative (sulfonamide charge/bond params).

### 5. guanidinium — TOO LOW (more negative, −37.4)
Driven by the guanidinium **charge** mis-distribution (central C charge +1.00 vs
RDKit +0.55) inflating electrostatics. **Fix path:** guanidinium formal-charge
rule.

## Charge outliers (mean |Δ| = 0.0092, but 23 atoms > 0.15)
- **ammonium** N: WebMM −1.74 vs RDKit −0.80 (Δ 0.94) — largest single charge error.
- **purine/guanine** ring N/C: Δ up to 0.47 — purine charge model.
- **guanidinium** central C: +1.00 vs +0.55.
- **carbon_disulfide** S: −0.28 vs −0.63.

## Typing gaps (51 atoms, 3.85%) — separate from energy, see CODE_STATUS
C5A/C5B alpha-beta (13), H_NAM/H_N3 (16), HNR+ charged-amine H (5),
cyclopropane CR3R (3), water OH2-vs-O_3 (3), N_2-vs-N_PL3 (3), CS₂ S2CM (2),
H_COOH-vs-H_ONC (2), minor (4).

## Summary
Typing and charges are ~96% / ~99% correct; the residual energy gap is dominated
by **5-ring angle parameters** (8 of 15 worst outliers) and **strained-ring
angles**, plus sparse P/S params and a few charge bugs. These are parameter-
table completeness issues, not algorithmic.
