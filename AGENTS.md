# Repository Agent Rules

This AGENTS.md extends the GLOBAL AGENTS.md.
All rules in this file are mandatory for this repository and may only
add stricter constraints, never weaken global rules.

Failure to follow these rules is an error.

---

## Stack Overview

Core library:
- Rust (edition 2021), `cargo` as the sole build/dependency manager
- `wasm-bindgen`/`wasm-pack` for the WASM target (demo entry: `site/index.html`,
  staged into the gitignored `pkg/`)
- MMFF94/MMFF94s force field, ETKDG v3 embedding, L-BFGS optimizer, MD engine,
  metadynamics, GBSA implicit solvation

Validation tooling:
- Python 3 scripts under `scripts/` (RDKit reference generation, benchmark/audit)
- RDKit is used only by scripts; the library itself has no Python runtime dependency

---

## Mandatory Workflow (In Addition to Global Rules)

For every non-trivial task:

1. Read `CODE_STATUS.md`
2. Create or overwrite `PLAN.md`
3. Implement strictly according to `PLAN.md`
4. Update `CODE_STATUS.md` after completion
5. CODE_STATUS.md and PLAN.md must follow their templates exactly. Do not add sections or free-form commentary.

Do not skip steps.
If a step cannot be completed, STOP and ask.

---

## Rust Dependency Management (Strict)

This repository uses **cargo** as the sole Rust dependency and build manager.

### Allowed Commands
- Add dependencies: `cargo add <crate>`
- Remove dependencies: `cargo remove <crate>`
- Build/test: `cargo build` / `cargo test`
- Lint: `cargo clippy --all-targets`
- WASM build: `wasm-pack build --target web --release`, then stage
  `site/index.html` into `pkg/` as the landing page

### Forbidden Commands
- `cargo install` for project-local tooling (use the pinned wasm-pack devDependency)
- Manually editing `Cargo.toml`/`Cargo.lock` unless explicitly instructed

### Rationale
- `cargo add/remove` are the canonical, state-aware workflow
- `Cargo.lock` is committed and must stay consistent with `Cargo.toml`

---

## Running Code & Tests

Examples:
- Run tests: `cargo test`
- Run lints: `cargo clippy --all-targets` (must be 0 warnings)
- Run MMFF benchmark: `python3 scripts/benchmark_mmff.py --no-speed`
- Run a diagnostic: `cargo run --release --example <name>`
- Build WASM: `wasm-pack build --target web --release && cp site/index.html pkg/index.html`

---

## Rust Coding Rules

- Keep `cargo clippy --all-targets` at 0 warnings; no dead code (use
  `#[allow(dead_code)]` only with a comment explaining why)
- No `unsafe` blocks
- Keep tests green: `cargo test` must pass with no flakes (pin RNG seeds in
  geometry/embedding tests — the ETKDG default seed is intentionally random)
- Follow existing module/error-handling patterns; no speculative refactors

---

## WASM & Site Rules

- `site/index.html` is the committed demo; `pkg/` is gitignored build output
- WASM exports in `src/lib.rs` are a public contract — do not break them
- New engine features that need browser access must be explicitly wired into
  `src/lib.rs` with tests, and reflected in `pkg/webmm.d.ts` after a build

---

## API Contract Rules

- WASM export signatures must remain stable
- Native changes must respect existing API contracts
- Breaking API changes require:
  - Explicit mention in `PLAN.md`
  - Documentation in `CODE_STATUS.md`

---

## Scope & Safety Constraints

- No refactors unless explicitly requested
- No dependency changes without approval
- No UI redesigns unless requested
- No speculative optimizations
- No “while we’re here” changes

---

## Completion Checklist

A task is complete only if:
- Implementation matches `PLAN.md`
- `CODE_STATUS.md` is updated accurately
- Changes are summarized clearly
