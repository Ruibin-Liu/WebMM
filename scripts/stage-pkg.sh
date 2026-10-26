#!/usr/bin/env bash
# stage-pkg.sh — build the WASM and stage the full deployable site into pkg/.
#
# Mirrors what .github/workflows/pages.yml stages in CI so a local pkg/
# preview matches the deployed Pages layout. wasm-pack never cleans its
# out-dir, so stale files from earlier eras (spike.html, old CDN-referencing
# pages...) accumulate otherwise — this prunes everything first.
#
# Staged layout (deploy root = pkg/):
#   index.html, playground.html   <- site/
#   vendor/3Dmol-min.js           <- app/vendor/   (site pages load vendored)
#   app/                          <- app/          (workbench + platform)
#   webmm.js, webmm_bg.wasm, ...  <- wasm-pack output
set -euo pipefail
cd "$(dirname "$0")/.."

if [ ! -d pkg ]; then
  echo "pkg/ does not exist; creating (wasm-pack will fill it)"
  mkdir -p pkg
fi

# prune everything wasm-pack will not rewrite (it never cleans its out-dir)
find pkg -mindepth 1 -maxdepth 1 -exec rm -rf {} +

npx wasm-pack build --target web --release

cp site/index.html pkg/index.html
cp site/playground.html pkg/playground.html
mkdir -p pkg/vendor
cp app/vendor/3Dmol-min.js pkg/vendor/3Dmol-min.js
cp -R app pkg/app

echo "staged: index, playground, vendor/3Dmol, app/ -> pkg/ ($(du -sh pkg | cut -f1))"
