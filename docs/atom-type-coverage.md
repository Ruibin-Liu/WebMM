# MMFF Atom Type Coverage

## Summary

| Metric | Count |
|--------|-------|
| Total MMFF types | 99 |
| Real types (crd > 0) | 86 |
| Pseudo/LP types (crd = 0) | 13 |
| **Currently covered** | **66 / 86 real types (76.7%)** |
| **Candidates found (ready to add)** | **7 more → 73/86 (84.9%)** |
| Remaining unfound | 9 real types |
| Molecules validated | 179 (all match RDKit < 0.01 kcal/mol) |

## Validation Sets

| Set | Molecules | Location |
|-----|-----------|----------|
| Original | 130 | `scripts/val_set/` |
| New (diverse chemistry) | 41 | `scripts/val_set_new/` |
| New (exotic types) | 8 | `scripts/val_set_new2/` |

## 66 Covered Atom Types

```
1  2  3  4  5  6  7  8  9     11      15 16 17 18 19 20
21 22 23 24 25 26 27 28       32      34 35 36 37 38    40
                       43 44 46 47            53      55 57
59 60 61      63 64 65            71      73 74 76 77 78
```

## 7 Types with Found Candidates (ready to add to validation set)

| Type | Element | SMILES | Molecule | Description |
|------|---------|--------|----------|-------------|
| 31 | H | `O` | Water | H on water O |
| 51 | O | `[o+]1ccc1` | Furanium cation | Aromatic oxonium O in 5-ring |
| 58 | N | `[n+]1(C)ccccc1` | N-methylpyridinium | N-methylated aromatic N |
| 68 | N | `[N+](C)(C)(C)[O-]` | Trimethylamine N-oxide | Quaternary N bonded to O⁻ |
| 69 | N | `O=n1ccccc1` | Pyridine N-oxide | Aromatic N-oxide N |
| 70 | O | `O` | Water | Water oxygen |
| 75 | P | `c1cc[p]1` | Aromatic phosphirene | P in aromatic 3-membered ring |

Found via RDKit test data (`bulk.sdf`) and systematic SMILES search.
Types 31+70 are already handled in WebMM (OH2/H_OH types) but not yet
exercised by validation molecules.

## 9 Remaining Unfound Types

| Type | Elem | Crd | Val | mltb | arom | sbmb | Description | Search notes |
|------|------|-----|-----|------|------|------|-------------|-------------|
| 48 | N | 2 | 2 | 0 | 0 | 0 | Divalent N (=N-, no H) | Azo/diazo/hydrazine all gave 9/47/8 |
| 52 | H | 1 | 1 | 0 | 0 | 0 | H variant | Water H gave 31, hydroxide untested |
| 54 | N | 3 | 4 | 2 | 0 | 1 | N in 3-ring with C=N | Diazirine variants gave 9/10 |
| 56 | N | 3 | 3* | 1 | 0 | 0 | N with mltb=1 | Enamine/amidine/hydroxylamine gave 8/40 |
| 67 | N | 3 | 4 | 2 | 0 | 1 | N in 4-ring with C=N | Azetine gave 20/30/40 |
| 79 | N | 2 | 3 | 2 | 1 | 0 | Aromatic N in diazine | Pyrazine/pyrimidine/phthalazine all gave 38 |
| 80 | C | 3 | 4 | 2 | 0 | 1 | C in 4-ring with C=C | Cyclobutene gave 30 (CE4R) |
| 81 | N | 3 | 4 | 1 | 1 | 1 | Aromatic N in small ring | Triazole/tetrazole gave 39/65/66 |
| 82 | N | 3 | 4 | 1 | 1 | 0 | Aromatic N variant | Pyrrole/indole/indazole gave 39 |

\* Type 56 has val=34 in MMFFProp, likely encoding artifact.

These 9 types require deeper investigation of RDKit's atom typing rules.
The sbmb types (54, 67, 80, 81) all involve small rings with multiple bonds.
Types 79 and 82 may require specific fused-ring substitution patterns.

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).

## Key Findings

### Types found via RDKit test data (bulk.sdf) + SMILES search (7 types)
- **Type 31** (H on water) + **Type 70** (water O): Simply add H₂O to validation set
- **Type 51** (aromatic oxonium O): 5-membered ring with O⁺ (furanium)
- **Type 58** (N-methylpyridinium N): N-alkylated aromatic ring
- **Type 68** (quaternary N-O): Trimethylamine N-oxide
- **Type 69** (pyridine N-oxide N): Aromatic N-oxide
- **Type 75** (aromatic phosphirene P): P in aromatic 3-ring

### Types that remain unfound (9 types)
These require deeper investigation of RDKit's internal atom typing algorithm.
The sbmb flag types (54, 67, 80, 81) involve exotic small-ring chemistry.
Types 79/82 may require specific fused-ring substitution patterns not yet tried.
