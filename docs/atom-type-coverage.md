# MMFF Atom Type Coverage

## Summary

| Metric | Count |
|--------|-------|
| Total MMFF types | 99 |
| Real types (crd > 0) | 86 |
| Pseudo/LP types (crd = 0) | 13 |
| **Currently covered** | **73 / 86 real types (84.9%)** |
| **Candidates found (ready to add)** | **7 more → 80/86 (93.0%)** |
| Remaining unfound | 2 real types |
| Molecules validated | 185 (all match RDKit < 0.01 kcal/mol) |

## Validation Sets

| Set | Molecules | Location |
|-----|-----------|----------|
| Original | 130 | `scripts/val_set/` |
| New (diverse chemistry) | 41 | `scripts/val_set_new/` |
| New (exotic batch 1) | 8 | `scripts/val_set_new2/` |
| New (exotic batch 2) | 6 | `scripts/val_set_new3/` |

## 7 Types with Found Candidates (ready to add)

| Type | Element | SMILES | Molecule | Source |
|------|---------|--------|----------|--------|
| 48 | N | `S(=N)(=O)C` | Sulfinylamine | RDKit source: "Divalent N replacing O in SO2" |
| 54 | N | `C=[NH2+]` | Methaniminium | RDKit source: "Iminium nitrogen N+=C" |
| 56 | N | `CN(C)C(=[NH2+])N(C)C` | Tetramethylguanidinium | RDKit source: "Guanidinium nitrogen" |
| 67 | N | `O=[n+]1cccc1` | Pyridine N-oxide (sp2 variant) | RDKit source: "sp2 N-oxide nitrogen" |
| 80 | C | `CN1C=C[N+](C)=C1` | N,N'-Dimethylimidazolium | RDKit source: "C between N's in imidazolium" |
| 81 | N | `CN1C=C[N+](C)=C1` | N,N'-Dimethylimidazolium | RDKit source: "Positive N in 5-ring" |
| 82 | N | `O=n1ccoc1` | Isoxazole N-oxide | RDKit source: "N-oxide in 5-ring" |

## 2 Remaining Unfound Types

| Type | Elem | Crd | Val | mltb | arom | sbmb | RDKit source comment | Search notes |
|------|------|-----|-----|------|------|------|---------------------|--------------|
| 52 | H | 1 | 1 | 0 | 0 | 0 | "HO=+ — H on oxenium oxygen" | Requires O+ type 51 with bonded H in ring. Extremely rare/unstable species. |
| 79 | N | 2 | 3 | 2 | 1 | 0 | "N5 — General N in 5-ring with alpha+beta heteroatoms" | Requires specific multi-heteroatom 5-ring where N has heteroatoms at both adjacent and non-adjacent positions. Existing NPYL/N5A/N5B classification may mask this type. |

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).

## Methodology

Types were found by:
1. Parsing RDKit's MMFFProp table for all 99 type definitions
2. Reading RDKit's `AtomTyper.cpp` source code for assignment conditions
3. Systematic SMILES search with RDKit MMFF typing
4. Cross-referencing with RDKit's `bulk.sdf` test data
5. Using `GetMMFFBondStretchParams` / `GetMMFFStretchBendParams` APIs for exact params
