# CDP regression suites

Headless-Chromium (Playwright) end-to-end suites against the two app pages.
All of them need:

- the test server on http://localhost:8901 serving the repo root
  (`python3 -m http.server 8901` from the repo root — `tests/cdp/run.sh`
  starts one if missing and only kills what it started)
- a Chromium binary at the Playwright path hard-coded in each file
  (update after Playwright cache upgrades)

Or just run everything, in dependency order:

    tests/cdp/run.sh            # all suites (regenerates fixtures if missing)
    tests/cdp/run.sh m6_platform m5_search_e2e   # a subset

## Suites

| file | page | covers |
|---|---|---|
| `m0_core` | index | 2D parse/render, properties, buttons, history, JSME |
| `m1_3d_pipeline` | index | embed → optimize → GFN-FF → exports (incl. demo parity) |
| `m2_conformers` | index | ensemble ranking, chart, ΔE reopt, cancel path |
| `m3_site_nav` | index | URL deep-links, history modal, axe a11y, keyboard nav |
| `m4_batch` | index | batch descriptors/QED/PAINS parity vs batch_refs.json, sort/filter/error rows, CSV content, 3D + SDF provenance |
| `m5_search_e2e` | index | search tab: library persistence, similarity/substructure parity, shape+pharm filter, RGD, scaffold frequency, pharmacophore screen |
| `m6_platform` | platform | LBDD platform page: shell/overview, rounds DAG + hit-as-query, consensus + provenance CSV, project export/import + multi-tab, triage overlay, neighborhood explorer funnel, accessibility (axe/keyboard/modal focus) |
| `m6_import_loop` | platform | import→reload race stress ×LOOP_N (snapshot guard + projection invalidation) |
| `m6_restore_rdkit_race` | platform | search-library restore racing the RDKit WASM load (deterministic: holds the wasm fetch) |
| `m6_csv_export` | both | CSV formula-injection neutralization + RFC 4180 quoting (unit cases per page + a real download flow) |
| `m6_sa_warmup` | platform | SA-table idle warmup firing before the engine is ready (holds the engine script; red/green) |
| `m6_analog_cancel` | platform | analog worker lifecycle: engine-fetch failure fails fast, superseded runs complete |
| `m6_analog_pool` | platform | persistent M3 worker pool: reuse across runs, poison recycling + job re-dispatch, pool health after chaos |
| `m6_analog_ui` | platform | explorer UI: thumbnail rendering at display scale, process-label tooltips, results CSV/SDF export |
| `m6_projects` | platform | M4 P1 project management: registry, per-project DBs, zero-copy migration, URL active project, export v1 envelope, import-as-new |
| `m6_analog_history` | platform | M4 P2 explorer persistence: auto-commit with lineage, history restore across reload, hide-as-undo |
| `m6_repoint` | platform | M4 P2 §3 pointer repoint: restore-into-current swaps the DB atomically, other tabs get the forced-reload modal, old DB trashed |
| `m6_scratch` | platform | M4 P3 scratch pad: legacy-history migration to webmm-scratch, IDB-backed clipboard, send-to-library gesture |

Note: since the workbench swap, m0–m5 run against `app/index.html` (the
former platform page, classic UI, full LBDD surface); m6* run against
`app/platform.html`.

## Fixture dependencies

Only two real inputs; everything else under /tmp is written by the suites
(downloads, screenshots) and can be deleted freely.

- **`/tmp/caff24.sdf`** — m1 only. 24-atom caffeine, regenerated via RDKit:
  `AddHs('Cn1cnc2c1c(=O)n(C)c(=O)n2C')` → ETKDG(42) → MMFF94s opt.
  Vanishes on reboot; `run.sh` regenerates it when missing (needs `rdkit`
  importable by python3).
- **`/tmp/m5-project.json`** — exported by a full `m6_platform` run and
  consumed by `m6_import_loop` and `m6_restore_rdkit_race`. Run m6 first
  (run.sh's default order does this).

## Environment variables

- `LOOP_N` (m6_import_loop) — iterations, default 10 (use 50 to hunt flakes)
- `STOP_ON_FAIL` (m6_import_loop) — stop at first failure with store/projection
  diagnostics dumped
- `KEEP=1` (scripts/build-ketcher.sh, not a test) — keep the build workspace

## Related

- `node tests/platform/run.js` — Node-only platform store suites (no browser)
- `scripts/build-ketcher.sh` — reproducible Ketcher bundle build (pinned
  versions; `--install` writes into app/vendor/ketcher)
