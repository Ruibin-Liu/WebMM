#!/usr/bin/env bash
# build-ketcher.sh — reproducible build of app/vendor/ketcher/ketcher-bundle.js.
#
# The committed bundle was a one-off esbuild IIFE of the official npm packages
# (no CDN at runtime; Indigo worker+wasm are embedded by ketcher-standalone's
# own dist). This script pins every input so the artifact can be rebuilt or
# version-bumped on demand:
#
#   ketcher-react      3.18.0   (React UI; Editor component)
#   ketcher-standalone 3.18.0   (Indigo WASM service, worker inlined; indigo 1.46.x)
#   esbuild            0.25.10  (bundler)
#
# Usage:
#   scripts/build-ketcher.sh              # build to a temp dir, print path
#   scripts/build-ketcher.sh --install    # also copy into app/vendor/ketcher/
#                                         # (ketcher-bundle.js, ketcher.css, LICENSE)
#
# Notes:
# - The banner guard only installs a `process` shim when the page has none —
#   platform.html defines its own global process(), and an unconditional
#   banner would silently overwrite it (documented pitfall).
# - ketcher.css is ketcher-react's dist stylesheet, copied unmodified.
# - Needs network for `npm install` (or a warm npm cache); everything else is local.
set -euo pipefail

PIN_KETCHER_REACT=ketcher-react@3.18.0
PIN_KETCHER_STANDALONE=ketcher-standalone@3.18.0
PIN_ESBUILD=esbuild@0.25.10
# browser shims for node builtins pulled in by ketcher-core (events/buffer/…)
PIN_SHIMS="events@3.3.0 buffer@6.0.3 process@0.11.10 util@0.12.5 path-browserify@1.0.1 stream-browserify@3.0.0"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENDOR="$ROOT/app/vendor/ketcher"
ENTRY="$ROOT/scripts/ketcher-entry.js"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ketcher-build.XXXXXX")"
# KEEP=1 preserves the workspace for inspection (smoke tests, diffing)
[ "${KEEP:-0}" = "1" ] || trap 'rm -rf "$WORK"' EXIT

echo "== workspace: $WORK"
cd "$WORK"
cat > package.json <<'EOF'
{ "name": "ketcher-build", "private": true }
EOF
# esbuild resolves imports from the ENTRY's directory — bundle a copy inside
# the workspace so $WORK/node_modules is found.
cp "$ENTRY" ./ketcher-entry.js
npm install --no-audit --no-fund --loglevel=error \
  $PIN_KETCHER_REACT $PIN_KETCHER_STANDALONE $PIN_ESBUILD $PIN_SHIMS

echo "== bundling (iife, minify)"
./node_modules/.bin/esbuild ./ketcher-entry.js \
  --bundle --format=iife --minify --platform=browser \
  --alias:path=path-browserify --alias:stream=stream-browserify \
  --define:process.env.NODE_ENV='"production"' \
  --banner:js=";(function(){if(typeof process==='undefined'){window.process={env:{NODE_ENV:'production'}}}else if(!process.env){process.env={NODE_ENV:'production'}}})();var global=globalThis;" \
  --outfile="$WORK/ketcher-bundle.js"

cp ./node_modules/ketcher-react/dist/index.css "$WORK/ketcher.css"
# neither npm package ships a LICENSE file at its root; vendor the Apache-2.0
# text once — reuse the one already committed when installing.
echo "== built artifacts"
ls -lh "$WORK"/ketcher-bundle.js "$WORK"/ketcher.css

echo "== sanity: banner guard + WebMMKetcher API present"
head -c 200 "$WORK/ketcher-bundle.js" | grep -q "typeof process==='undefined'" || { echo "FAIL: banner guard missing"; exit 1; }
grep -q 'WebMMKetcher' "$WORK/ketcher-bundle.js" || { echo "FAIL: WebMMKetcher missing"; exit 1; }
for api in mount unmount getMolfile setMolecule; do
  grep -q "$api(" "$WORK/ketcher-bundle.js" || { echo "FAIL: API $api missing"; exit 1; }
done
echo "OK"

if [ "${1:-}" = "--install" ]; then
  cp "$WORK/ketcher-bundle.js" "$VENDOR/ketcher-bundle.js"
  cp "$WORK/ketcher.css"       "$VENDOR/ketcher.css"
  # LICENSE: packages don't ship one; keep the committed Apache-2.0 text.
  echo "== installed into $VENDOR (re-run the page Ketcher smoke test after)"
else
  echo "== dry run only; pass --install to write into app/vendor/ketcher/"
fi
