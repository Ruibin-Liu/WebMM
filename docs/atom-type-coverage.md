# MMFF Atom Type Coverage

## Summary

| Metric | Count |
|--------|-------|
| Total MMFF types | 99 |
| Real types (crd > 0) | 86 |
| Pseudo/LP types (crd = 0) | 13 |
| **Covered** | **66 / 86 real types (76.7%)** |
| Missing | 16 real types |
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

## 16 Uncovered Real Types

| Type | Elem | Crd | Val | mltb | arom | lin | sbmb | Description | Best Candidate | Status |
|------|------|-----|-----|------|------|-----|------|-------------|-----------------|--------|
| 31 | H | 1 | 1 | 0 | 0 | 0 | 0 | H variant (specific C context) | Unknown — allene, acetylene gave type 5 | ❌ Hard |
| 48 | N | 2 | 2 | 0 | 0 | 0 | 0 | Divalent N (=N- terminal, no H) | Unknown — azo/diazo gave 9 or 47 | ❌ Hard |
| 51 | O | 2 | 3 | 2 | 0 | 0 | 0 | O in cumulated system (O=C=C) | Unknown — ketene/CO2 gave type 7 | ❌ Hard |
| 52 | H | 1 | 1 | 0 | 0 | 0 | 0 | H variant (charged species?) | Unknown | ❌ Unknown |
| 54 | N | 3 | 4 | 2 | 0 | 0 | 1 | N in 3-ring with C=N | Unknown — diazirine gave 9/10 | ❌ Very hard |
| 56 | N | 3 | 3* | 1 | 0 | 0 | 0 | N with mltb=1 (val field unusual) | Unknown — enamine/amidine gave 8/40 | ❌ Unknown |
| 58 | N | 3 | 4 | 1 | 1 | 0 | 1 | Aromatic N in N-methyl pyridinium | `[n+]1(C)ccccc1` | ✅ Found |
| 67 | N | 3 | 4 | 2 | 0 | 0 | 1 | N in 4-ring with C=N | Unknown — azetine gave 20/30/40 | ❌ Very hard |
| 68 | N | 4 | 4 | 0 | 0 | 0 | 0 | Quaternary ammonium (specific) | Unknown — NMe4+ gave type 34 | ❌ Hard |
| 69 | N | 3 | 4 | 1 | 1 | 0 | 0 | Aromatic N (pyridinium-like) | Unknown — pyridinium gave 38 | ❌ Hard |
| 70 | O | 2 | 2 | 0 | 0 | 0 | 0 | Ether O variant | Unknown — all ethers gave type 6 | ❌ Hard |
| 75 | P | 2 | 3 | 2 | 0 | 0 | 1 | P in aromatic 3-ring | `c1cc[p]1` (aromatic phosphirene) | ✅ Found |
| 79 | N | 2 | 3 | 2 | 1 | 0 | 0 | Aromatic N in diazine | Unknown — pyrazine/quinoline gave 38 | ❌ Hard |
| 80 | C | 3 | 4 | 2 | 0 | 0 | 1 | C in 4-ring with C=C (variant) | Unknown — cyclobutene gave type 30 | ❌ Hard |
| 81 | N | 3 | 4 | 1 | 1 | 0 | 1 | Aromatic N in small ring | Unknown — triazole gave 39/65/66 | ❌ Very hard |
| 82 | N | 3 | 4 | 1 | 1 | 0 | 0 | Aromatic N (pyrrole-like variant) | Unknown — pyrrole/indazole gave 39 | ❌ Hard |

\* Type 56 has val=34 in the MMFFProp table, likely a packed/encoding artifact.

## 13 Pseudo/LP Types (not real atoms — untestable)

Types 87–99: lone-pair, dummy, and metal-ion pseudo-atoms (crd=0, val=0).

## Key Findings

### Types with found candidates (2)
- **Type 58**: N-methylpyridinium `[n+]1(C)ccccc1` — aromatic N in a methylated pyridine
- **Type 75**: Aromatic phosphirene `c1cc[p]1` — P in aromatic 3-membered ring

### Pattern: RDKit assigns more common types than expected (11 types)
Most missing types appear to be "specialized variants" that RDKit only assigns in
narrow contexts that our candidate molecules didn't trigger. For example:
- Type 68 (quaternary N): All ammonium compounds tested give type 34 (N_4)
- Type 69/79/82 (aromatic N): All aromatic N compounds give types 38 (N_AR) or 39 (NPYL)
- Type 70 (ether O): All ethers give type 6 (O_R)
- Type 80 (4-ring C): Cyclobutene gives type 30 (CE4R)

These may require very specific substitution patterns or ring fusions to trigger.

### Types that are genuinely exotic (3 types)
- **Type 54** (3-ring N with C=N): Requires a diazirine-like structure with specific bond representation
- **Type 67** (4-ring N with C=N): Requires a specific azetine isomer
- **Type 81** (aromatic N in 3-ring): Extremely exotic aromatic small-ring nitrogen

## Notes for Future Work

1. **Types 58 and 75** can be added immediately — candidate molecules identified
2. **Types 31, 48, 51, 52, 56** — need deeper investigation of RDKit's atom typing rules
   to understand what structural features trigger these types
3. **Types 54, 67, 80, 81** — very exotic small-ring chemistry, may require
   specialized synthetic precursors
4. **Types 68, 69, 70, 79, 82** — RDKit may assign these only in specific fused-ring
   or highly-substituted contexts; need systematic SMARTS-pattern investigation
5. The **sbmb** flag (small-ring multiple bond) appears on types 54, 58, 67, 75, 80, 81
   — all involve small rings (3-4 membered) with double bond character
