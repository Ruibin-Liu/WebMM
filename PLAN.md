# Plan: Phases 6-9 — Parameter Loading, Test Fixes, Clippy, Docs

## Task 1: Embed and load MMFF parameters from JSON
- Files: src/utils/mod.rs, src/mmff/bond.rs, data/mmff94_sample_parameters.json

### Task 2: Fix ignored parser tests
- Files: src/molecule/parser.rs

### Task 3: Add end-to-end integration test
- Files: src/lib.rs

### Task 4: Fix all clippy warnings
- Files: src/mmff/mod.rs, src/mmff/bond.rs, src/mmff/angle.rs, src/mmff/torsion.rs, src/mmff/vdw.rs, src/molecule/graph.rs, src/molecule/parser.rs, src/molecule/mod.rs, src/etkdg/mod.rs, src/optimizer/mod.rs, src/lib.rs, Cargo.toml

### Task 5: Fix pkg/index.html API mismatches
- Files: pkg/index.html

### Task 6: Update CODE_STATUS.md
- Files: CODE_STATUS.md, PLAN.md
