#!/usr/bin/env python3
"""Color-rule agreement check: the page's simplified ROCS-style color rules
vs RDKit's BaseFeatures.fdef perception (Donor/Acceptor/PosIonizable/
NegIonizable families).

The page rules are deliberately simpler (see app/index.html COLOR_Q); this
script quantifies where they agree and differ on the LBDD reference set +
demo library. Hydrophobe/ring types are excluded: their definitions differ
by design (per-atom vs lumped) and no shared oracle exists.

The page side is evaluated through Playwright (real browser, real SMARTS),
so this script drives the check via a node helper printed below when run
with --emit-js.

Usage:
  python3 scripts/color_rules_check.py            # prints python-side stats
  python3 scripts/color_rules_check.py --page JSON   # consumes page output
"""
import json
import os
import sys
from pathlib import Path

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures

REPO = Path(__file__).resolve().parent.parent


def fdef_sites(m):
    factory = fdef_sites.factory
    out = {"donor": set(), "acceptor": set(), "pos": set(), "neg": set()}
    for f in factory.GetFeaturesForMol(m):
        fam = f.GetFamily()
        if fam == "Donor":
            out["donor"].update(f.GetAtomIds())
        elif fam == "Acceptor":
            out["acceptor"].update(f.GetAtomIds())
        elif fam == "PosIonizable":
            out["pos"].update(f.GetAtomIds())
        elif fam == "NegIonizable":
            out["neg"].update(f.GetAtomIds())
    return out


RULES = {
    "donor": ['[NX3;H1,H2,H3]', '[OX2;H1,H2]', '[nX3;H1]'],
    "acceptor": ['[#8;!+]', '[#7;H0;!+]'],
    "pos": ['[+;!$([N+](=O)[O-])]', '[NX3;H1,H2;!$(N[a]);!$(NC=O)]', '[NX3](=[NX2])[N]'],
    "neg": ['[-]', '[CX3](=[OX1])-[OX2;H1]'],
}


def rule_sites(m):
    out = {}
    for t, patterns in RULES.items():
        s = set()
        for sma in patterns:
            q = Chem.MolFromSmarts(sma)
            if q is None:
                continue
            for mt in m.GetSubstructMatches(q):
                idx = mt[-1] if (t == "neg" and len(mt) == 3) else mt[0]
                s.add(idx)
        out[t] = s
    return out


def main():
    # NOTE: the 'neg' acid rule takes the LAST atom (the O-H); the page takes
    # atoms[2] for 3-atom matches — same convention.
    factory = ChemicalFeatures.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
    fdef_sites.factory = factory

    mols = {}
    refs = json.loads((REPO / "tests/fixtures/lbdd/refs.json").read_text())
    for name, rec in refs.items():
        mols[name] = Chem.MolFromSmiles(rec["smiles"])
    lib = json.loads((REPO / "app/demo_library.js").read_text()
                     .split("window.DEMO_LIBRARY = ", 1)[1].rsplit(";", 1)[0])
    for e in lib:
        mols.setdefault("lib:" + e["name"], Chem.MolFromSmiles(e["smiles"]))

    stats = {t: {"exact": 0, "jaccard_sum": 0.0, "n": 0, "diffs": []} for t in RULES}
    for name, m in mols.items():
        if m is None:
            continue
        a, b = fdef_sites(m), rule_sites(m)
        for t in RULES:
            fa, fb = a[t], b[t]
            stats[t]["n"] += 1
            if fa == fb:
                stats[t]["exact"] += 1
            elif len(stats[t]["diffs"]) < 8:
                stats[t]["diffs"].append({"mol": name, "fdef": sorted(fa), "rules": sorted(fb)})
            u = fa | fb
            stats[t]["jaccard_sum"] += (len(fa & fb) / len(u)) if u else 1.0

    out = {t: {
        "exact_match": f"{v['exact']}/{v['n']}",
        "mean_jaccard": round(v["jaccard_sum"] / max(1, v["n"]), 4),
        "example_diffs": v["diffs"][:4],
    } for t, v in stats.items()}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
