#!/usr/bin/env python3
"""Generate tests/fixtures/lbdd/search_refs.json — Python-RDKit gold references
for the Workbench Search tab (similarity Tanimoto + substructure hits).

Queries cover both modes: plain-SMILES similarity queries and SMARTS
substructure queries. Library = the bundled demo library (same source as
app/demo_library.js, re-derived here for independence).

Reference semantics:
  - similarity: Tanimoto over the six parity-validated bit fingerprints
    (Morgan r2/2048, RDKit 2048, MACCS 167, AtomPair 2048, TopologicalTorsion
    2048), full float precision (JS double == C++ double on integer ratios)
  - substructure: Chem.HasSubstructMatch defaults (no AddHs) — the standard
    RDKit semantics the page replicates

Usage: python3 scripts/gen_search_refs.py
"""
import json
import sys
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "tests" / "fixtures" / "lbdd" / "search_refs.json"
LIB_SCRIPT = REPO / "app" / "demo_library.js"

QUERIES = {
    # name -> {"smiles": ..., "kind": "sim"|"sub"}   (sub queries may be SMARTS)
    "aspirin": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "kind": "sim"},
    "benzene": {"smiles": "c1ccccc1", "kind": "sim"},
    "pyridine_ring": {"smiles": "c1ccncc1", "kind": "sub"},
    "phenol": {"smiles": "Oc1ccccc1", "kind": "sub"},
    "carboxylic": {"smiles": "C(=O)O", "kind": "sub"},
    "amine_smarts": {"smiles": "[NX3;H2,H1;!$(NC[S,O]=O)]", "kind": "sub"},   # primary/secondary non-acylated amine
}


def fp_generators():
    from rdkit.Chem import MACCSkeys
    return {
        "morgan": lambda m: rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(m),
        "rdkit": lambda m: rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048).GetFingerprint(m),
        "maccs": MACCSkeys.GenMACCSKeys,
        "atompair": lambda m: rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048).GetFingerprint(m),
        "topologicaltorsion": lambda m: rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048).GetFingerprint(m),
    }


def load_library():
    text = LIB_SCRIPT.read_text()
    payload = text.split("window.DEMO_LIBRARY = ", 1)[1].rsplit(";", 1)[0]
    return json.loads(payload)


def main() -> int:
    lib = load_library()
    gens = fp_generators()
    mols, valid = [], []
    for e in lib:
        m = Chem.MolFromSmiles(e["smiles"])
        if m is None:
            continue
        mols.append(m)
        valid.append(e["name"])

    fps = {k: [] for k in gens}
    for key, gen in gens.items():
        fps[key] = [gen(m) for m in mols]

    refs = {"library": valid, "queries": {}}
    for qname, q in QUERIES.items():
        qmol = Chem.MolFromSmiles(q["smiles"])
        if qmol is None:
            qmol = Chem.MolFromSmarts(q["smiles"])
        assert qmol is not None, qname
        rec = {"kind": q["kind"], "smiles": q["smiles"]}
        if q["kind"] == "sim":
            rec["tanimoto"] = {
                key: list(DataStructs.BulkTanimotoSimilarity(gen(qmol), vec))
                for key, gen, vec in ((k, g, fps[k]) for k, g in gens.items())
            }
        else:
            rec["hits"] = sorted(n for n, m in zip(valid, mols) if m.HasSubstructMatch(qmol))
        refs["queries"][qname] = rec

    OUT.write_text(json.dumps(refs, indent=1))
    n_hits = {q: len(r.get("hits", [])) for q, r in refs["queries"].items()}
    print(f"wrote {OUT}: library={len(valid)}, queries={len(refs['queries'])}, hits={n_hits}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
