#!/usr/bin/env bash
# run.sh — run the CDP suites in dependency order.
#
#   tests/cdp/run.sh                 # everything (m0-m6, platform Node, 专项)
#   tests/cdp/run.sh m6_platform m5  # a subset (names are prefix-matched)
#
# Handles the two fixture dependencies automatically:
#   /tmp/caff24.sdf      regenerated via RDKit when missing (m1 input)
#   /tmp/m5-project.json produced by the full m6_platform run before the
#                        suites that consume it (default order only)
# A server is started on :8901 only if none is listening, and only a server
# started here is stopped here.
set -u
cd "$(dirname "$0")"

PORT=8901
WE_STARTED_SERVER=0

ensure_server() {
  if curl -s -o /dev/null --max-time 2 "http://localhost:$PORT/app/platform.html"; then
    echo "== server already up on :$PORT"
    return 0
  fi
  echo "== starting test server on :$PORT (from the repo root)"
  ( cd ../.. && python3 -m http.server $PORT >/tmp/webmm-cdp-server.log 2>&1 & )
  WE_STARTED_SERVER=1
  for _ in 1 2 3 4 5 6 7 8 9 10; do
    curl -s -o /dev/null --max-time 1 "http://localhost:$PORT/app/platform.html" && return 0
    sleep 0.5
  done
  echo "FAIL: test server did not come up (see /tmp/webmm-cdp-server.log)"; exit 1
}

ensure_caff24() {
  if [ -s /tmp/caff24.sdf ]; then return 0; fi
  echo "== regenerating /tmp/caff24.sdf via RDKit (reboot wiped it)"
  python3 - <<'PY' || { echo "FAIL: could not regenerate caff24.sdf (is rdkit importable?)"; exit 1; }
from rdkit import Chem
from rdkit.Chem import AllChem
m = Chem.AddHs(Chem.MolFromSmiles('Cn1cnc2c1c(=O)n(C)c(=O)n2C'))
AllChem.EmbedMolecule(m, randomSeed=42)
AllChem.MMFFOptimizeMolecule(m)
print(Chem.MolToMolBlock(m), file=open('/tmp/caff24.sdf', 'w'))
assert m.GetNumAtoms() == 24
PY
}

stop_own_server() {
  if [ "$WE_STARTED_SERVER" = "1" ]; then
    pkill -f "http.server $PORT" 2>/dev/null && echo "== stopped the server we started"
  fi
}
trap stop_own_server EXIT

ALL_ORDER="m0_core m1_3d_pipeline m2_conformers m3_site_nav m4_batch m5_search_e2e m6_platform m6_import_loop m6_restore_rdkit_race m6_csv_export m6_sa_warmup m6_analog_cancel m6_analog_pool m6_analog_ui"

if [ $# -gt 0 ]; then
  SUITES=""
  for want in "$@"; do
    hit=""
    for t in $ALL_ORDER; do
      case "$t" in
        "$want"*) hit="$t"; break ;;
        *_"$want"*) [ -z "$hit" ] && hit="$t" ;;   # m6 matches m6_* (first hit wins above)
      esac
    done
    if [ -z "$hit" ]; then echo "unknown suite: $want (candidates: $ALL_ORDER)"; exit 2; fi
    SUITES="$SUITES $hit"
  done
else
  SUITES="$ALL_ORDER"
fi

ensure_server
for t in $SUITES; do
  case "$t" in m1_*) ensure_caff24 ;; esac
done

pass=0; fail=0; failed=""
for t in $SUITES; do
  printf '%-22s ' "$t"
  if node "$t.test.js" > "/tmp/webmm-$t.log" 2>&1; then
    tail -1 "/tmp/webmm-$t.log"; pass=$((pass+1))
  else
    echo "FAIL (log: /tmp/webmm-$t.log)"; fail=$((fail+1)); failed="$failed $t"
  fi
done

echo "== node platform suite"
printf '%-22s ' "tests/platform/run.js"
if ( cd ../.. && node tests/platform/run.js > /tmp/webmm-platform-node.log 2>&1 ); then
  tail -1 /tmp/webmm-platform-node.log; pass=$((pass+1))
else
  echo "FAIL (log: /tmp/webmm-platform-node.log)"; fail=$((fail+1)); failed="$failed platform-node"
fi

echo "RESULT: $pass suites passed, $fail failed${failed:+ —failed:$failed}"
[ "$fail" -eq 0 ]
