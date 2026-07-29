# Plan: Fix discrepancies found in new validation compounds

## Final Status: 41/41 types match, 32/41 energy <0.01

## Completed
All atom typing patterns fixed (13 new type detections):
- Sulfonium/hypervalent S (17, 18), Cumulated N (53, 47), Guanidinium (55, 57, 36)
- H on P (71), Nitroso N (46), Isonitrile (60, 61), Nitrile oxide O (35)
- Amide/imine N guard, Nitro charge guard

Algorithm fixes:
- Torsion lookup: central atoms now use eq levels
- CR3R eq levels: [22,22,22,1] → [22,22,20,1]
- Guanidinium charge distribution (+1/n_nitrogens)

~40 bond/angle entries added from RDKit verbose.

## Remaining (9 molecules, all <0.75 kcal/mol)
All param-precision level (empirical rule 3-dp rounding):
- methyl_nitrite (+0.73): torsion V2 factor
- SF4 (-0.66): S-F empirical rule precision
- methyl_isothiocyanate (+0.61): C=S empirical rule
- methyl_azide (+0.39): N=N empirical rule
- methyl_isocyanate (+0.39): cumulated bond precision
- guanidinium_salt (+0.27): CGD+-NCN+ bond precision
- 3 more <0.15
