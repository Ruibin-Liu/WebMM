# Plan: GitHub CI/pages audit pass — rustfmt clean, pin wasm-pack in CI

## Context
Careful audit of the GitHub setup found:
1. **Rust CI is red on every recent commit** — `cargo fmt --all -- --check` fails
   with 505 hunks across 28 files (long-standing; clippy/test/build-wasm steps
   are skipped or green, but the lint job never passes). `cargo clippy -- -D
   warnings` itself passes locally.
2. **wasm-pack version drift CI vs local**: `pages.yml` and `rust.yml` install
   wasm-pack via the unpinned installer script (`curl init.sh | sh`) → CI gets
   v0.15.0, while local dev builds with the pinned devDependency 0.13.1
   (package.json `^0.13.0`). The deployed Pages artifact is therefore built
   with a different tool than local — reproducibility gap.
3. Pages deploys themselves are healthy: latest push (d48cfc1, incl. eb57769
   src changes) deployed successfully per the Actions API. Live page runtime
   cannot be verified from this environment (github.io is SSL-blocked here).

## Phase 1 — `cargo fmt --all` (formatting-only commit)
- Run `cargo fmt --all`; commit the result. Purely mechanical whitespace
  changes (same rustfmt 1.8.0/rustc 1.92.0 as CI's dtolnay@stable, so output
  matches). No semantics change; this unblocks the lint job's fmt step.

## Phase 2 — Pin wasm-pack v0.13.1 in CI workflows
- Replace the unpinned installer step in `.github/workflows/pages.yml` and
  `.github/workflows/rust.yml` with a pinned download of
  `wasm-pack-v0.13.1-x86_64-unknown-linux-musl.tar.gz` (verified to exist on
  the v0.13.1 release) → CI builds with the same version as local dev.

## Phase 3 — Docs per workflow
- Overwrite `PLAN.md` (this file).
- Update `CODE_STATUS.md`: Recently Completed entry + CI status note.

## Verification
- `cargo fmt --all -- --check` — 0 diffs.
- `cargo clippy -- -D warnings` — 0 warnings (exactly as CI's lint job runs).
- `cargo test` — 199/199.
- `wasm-pack build --target web --release` (local 0.13.1) — succeeds; staged
  `pkg/` demo is the same contract as CI output.
- Push; poll GitHub Actions until both workflows complete green.

## Out of scope
- Live-site runtime verification (github.io unreachable from this environment;
  user opens https://ruibin-liu.github.io/WebMM/ in a browser).
- Reformatting `scripts/` Python or `site/index.html` (not rustfmt's domain).
- Changing the deployed artifact layout or the pages.yml trigger paths.
