# Plan: Release v0.6.1 (patch)

## Context
Since v0.6.0 (922cfcb) the repo accumulated patch-level changes only: demo
page WASM API fixes (embed/optimize/metad buttons), Rust CI green (cargo fmt,
clippy 1.97 collapsible_match, wasm-pack v0.13.1 pin), README refresh, debug
print cleanup. No new features, no breaking API changes → semver patch bump
to 0.6.1, consistent with the repo convention of tagging at the crate
version.

## Phase 1 — Version bump
- `Cargo.toml`, `Cargo.lock` (webmm package), `package.json`: 0.6.0 → 0.6.1.
- Verify: `cargo check` (lock consistency), `cargo test` 199/199, clippy 0.

## Phase 2 — Docs per workflow
- Overwrite `PLAN.md` (this file).
- Update `CODE_STATUS.md` (Recently Completed).

## Phase 3 — Release
- Commit `release v0.6.1: ...`; tag `v0.6.1`; push commit + tag.
- Pages workflow fires on `v*` tags → demo redeploys (expected, same artifact).
- No GitHub Release created (repo has none; tags only).

## Out of scope
- Changelog file / release notes body (repo has no CHANGELOG; notes go in
  the commit message).
- Any feature or API changes (patch release only).
