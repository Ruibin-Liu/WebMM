#!/usr/bin/env python3
"""Generate app/lbdd_data.js — the client-side LBDD catalog for the Workbench.

Sources (both pinned & verified against the vendored RDKit minimal 2026.03.6):
  - PAINS A/B/C + BRENK structural-alert SMARTS:
      RDKit GitHub, tag Release_2026_03_6,
      Code/GraphMol/FilterCatalog/{pains_a,pains_b,pains_c,brenk}.in
      (C array literals: {"name", "smarts", 0, ""}; adjacent string literals
       are implicit concatenations and may span lines)
  - QED constants (Bickerton 2012 as implemented in RDKit):
      the pure-Python rdkit/Chem/QED.py shipped with the locally installed
      RDKit (weights x3, ADS parameters x8, AcceptorSmarts x11,
      StructuralAlertSmarts x117), parsed with the ast module.

Output is a committed artifact (offline-capable site); rerun to regenerate.

Usage:
  python3 scripts/gen_lbdd_data.py            # reads cache/urls as needed
  python3 scripts/gen_lbdd_data.py --refresh  # force re-download
"""
import argparse
import ast
import json
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "app" / "lbdd_data.js"

RDKIT_TAG = "Release_2026_03_6"
CATALOG_FILES = {  # output key -> .in file(s)
    "pains": ["pains_a.in", "pains_b.in", "pains_c.in"],
    "brenk": ["brenk.in"],
}
BASE_URL = (
    "https://raw.githubusercontent.com/rdkit/rdkit/"
    + RDKIT_TAG
    + "/Code/GraphMol/FilterCatalog/"
)

# The pure-Python QED implementation shipped with local RDKit.
QED_PY_CANDIDATES = [
    Path("/opt/homebrew/lib/python3.14/site-packages/rdkit/Chem/QED.py"),
]
try:  # fallback via import if the hardcoded path misses
    import rdkit

    QED_PY_CANDIDATES.append(Path(rdkit.__file__).parent / "Chem" / "QED.py")
except Exception:  # pragma: no cover
    pass

EXPECT = {
    "pains": 480,  # A=16 + B=55 + C=409 (matches FilterCatalogParams A+B+C)
    "brenk": 105,
    "qed_acceptors": 11,
    "qed_alerts": 116,
}

_STR = r'"(?:[^"\\]|\\.)*"'
ENTRY_RE = re.compile(
    (
        r'\{\s*(' + _STR + r'(?:\s*' + _STR + r')*)'
        r'\s*,\s*(' + _STR + r'(?:\s*' + _STR + r')*)'
        r'\s*,\s*\d+\s*,\s*""\s*\}'
    ),
    re.S,
)
_STR_PIECE_RE = re.compile(_STR)


def join_str_parts(group: str) -> str:
    """Concatenate the quoted pieces of one logical C string, unescaped."""
    return "".join(p[1:-1] for p in _STR_PIECE_RE.findall(group)).replace('\\"', '"')


def fetch_catalog(fname: str, refresh: bool) -> str:
    cache = Path(tempfile.gettempdir()) / f"webmm_lbdd_{RDKIT_TAG}_{fname}"
    if cache.exists() and not refresh:
        return cache.read_text()
    url = BASE_URL + fname
    print(f"fetching {url}")
    text = urllib.request.urlopen(url, timeout=60).read().decode()
    cache.write_text(text)
    return text


def parse_catalog(text: str) -> list:
    entries = []
    for m in ENTRY_RE.finditer(text):
        name, smarts = join_str_parts(m.group(1)), join_str_parts(m.group(2))
        if not smarts:  # placeholder rows in some .in files
            continue
        entries.append({"id": name, "smarts": smarts})
    return entries


def find_qed_py() -> Path:
    for cand in QED_PY_CANDIDATES:
        if cand and cand.exists():
            return cand
    raise SystemExit("rdkit/Chem/QED.py not found; install RDKit python package")


def parse_qed_py(path: Path) -> dict:
    tree = ast.parse(path.read_text())
    weights, ads, acceptors, alerts = {}, {}, None, None

    def const_str(node):
        """Evaluate a Constant or BinOp(Add)-concatenated string."""
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return node.value
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return const_str(node.left) + const_str(node.right)
        return ast.literal_eval(node)

    def const_list(node):
        return [const_str(e) for e in node.elts]

    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        name, value = target.id, node.value
        if name in ("WEIGHT_MAX", "WEIGHT_MEAN", "WEIGHT_NONE") and isinstance(value, ast.Call):
            weights[name.removeprefix("WEIGHT_").lower()] = [
                ast.literal_eval(a) for a in value.args
            ]
        elif name == "AcceptorSmarts" and isinstance(value, ast.List):
            acceptors = const_list(value)
        elif name == "StructuralAlertSmarts" and isinstance(value, ast.List):
            alerts = const_list(value)
        elif name == "adsParameters" and isinstance(value, ast.Dict):
            for k, v in zip(value.keys, value.values):
                assert isinstance(v, ast.Call)
                # ADSparameter(A=..., B=..., ...) — all keyword args
                kw = {a.arg: ast.literal_eval(a.value) for a in v.keywords}
                if not kw and v.args:  # tolerate positional form
                    kw = dict(zip(["A", "B", "C", "D", "E", "F", "DMAX"],
                                  [ast.literal_eval(a) for a in v.args]))
                ads[k.value] = kw
    if not weights or not ads or not acceptors or not alerts:
        raise SystemExit(f"QED.py parse incomplete: {path}")
    return {"weights": weights, "ads": ads, "acceptors": acceptors, "alerts": alerts}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()

    data = {}
    for key, files in CATALOG_FILES.items():
        entries = []
        for f in files:
            entries.extend(parse_catalog(fetch_catalog(f, args.refresh)))
        data[key] = entries

    qed_py = find_qed_py()
    data["qed"] = parse_qed_py(qed_py)
    # descriptor order for weights: MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS
    data["qed"]["ads_order"] = ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]

    # ---- self-checks (hard gate) ----
    assert len(data["pains"]) == EXPECT["pains"], f"pains {len(data['pains'])} != {EXPECT['pains']}"
    assert len(data["brenk"]) == EXPECT["brenk"], f"brenk {len(data['brenk'])} != {EXPECT['brenk']}"
    assert len(data["qed"]["acceptors"]) == EXPECT["qed_acceptors"]
    assert len(data["qed"]["alerts"]) == EXPECT["qed_alerts"]
    for w in data["qed"]["weights"].values():
        assert len(w) == 8
    assert set(data["qed"]["ads"]) == set(data["qed"]["ads_order"])
    for cat in ("pains", "brenk"):
        for e in data[cat]:
            assert e["id"] and e["smarts"], e

    provenance = {
        "catalog_source": f"rdkit {RDKIT_TAG} Code/GraphMol/FilterCatalog/*.in",
        "qed_source": str(qed_py),
        "counts": {
            "pains": len(data["pains"]),
            "brenk": len(data["brenk"]),
            "qed_alerts": len(data["qed"]["alerts"]),
            "qed_acceptors": len(data["qed"]["acceptors"]),
        },
    }

    js = (
        "/* Generated by scripts/gen_lbdd_data.py — do not edit by hand.\n"
        " * Provenance: "
        + json.dumps(provenance, indent=1).replace("\n", "\n * ")
        + "\n"
        " */\n"
        "window.LBDD_DATA = "
        + json.dumps(data, separators=(",", ":"))
        + ";\n"
    )
    OUT.write_text(js)
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")
    print(json.dumps(provenance["counts"]))


if __name__ == "__main__":
    sys.exit(main())
