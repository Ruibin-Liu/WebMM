# CDP regression suites

Headless-Chromium end-to-end suites (Playwright). They need:
- the test server on http://localhost:8901 serving the repo root
  (`python3 -m http.server 8901` from the repo root)
- a Chromium binary at the PLAYWRIGHT path hard-coded in each file
  (update after Playwright cache upgrades)
- /tmp/caff24.sdf for m1 (24-atom caffeine, regenerated via RDKit:
  AddHs('Cn1cnc2c1c(=O)n(C)c(=O)n2C') → ETKDG(42) → MMFF94s opt
  (24-atom caffeine; an earlier recipe here was theophylline (21 atoms))

Run all (from this directory):
    node m0_core.test.js && node m1_3d_pipeline.test.js && \
    node m2_conformers.test.js && node m3_site_nav.test.js && \
    node m4_batch.test.js && node m5_search_e2e.test.js
- m6_platform.test.js — the LBDD platform page (app/platform.html): rounds DAG + hit-as-query, consensus + provenance CSV, project export/import + multi-tab broadcast, triage overlay, neighborhood explorer (full funnel)

Suites:
- m0_core           2D parse/render, properties, buttons, history, JSME
- m1_3d_pipeline    embed → optimize → GFN-FF → exports (incl. demo parity)
- m2_conformers     ensemble ranking, chart, ΔE reopt, cancel path
- m3_site_nav       URL deep-links, history modal, axe a11y, keyboard nav
- m4_batch          batch descriptors/QED/PAINS parity vs batch_refs.json,
                    sort/filter/error-row/SDF-input/cancel/row-reload,
                    CSV content, 3D + SDF provenance
- m5_search_e2e     Search tab: library persistence, similarity/substructure
                    parity vs tests/fixtures/lbdd (exact), shape+pharm filter,
                    RGD, scaffold frequency, pharmacophore screen
