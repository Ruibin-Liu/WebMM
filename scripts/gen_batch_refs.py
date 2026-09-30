#!/usr/bin/env python3
"""Generate tests/fixtures/lbdd/batch_refs.json — Python-RDKit gold references
for the Workbench BATCH table (m4 CDP suite).

refs.json (LBDD) stores QED.properties — QED-flavored HBA/HBD counts that
differ from the Lipinski descriptors the batch table renders (e.g. aspirin:
QED HBA 4 vs NumHBA 3). This fixture therefore records the exact descriptor
set describeBatchMolecule displays:

  MW / cLogP / TPSA     average MW, Crippen MolLogP, TPSA
  HBD / HBA / RotB      rdMolDescriptors CalcNumHBD/CalcNumHBA/
                        CalcNumRotatableBonds
  qed_mean / qed_max    QED scores (6 dp)
  pains / brenk         matched-catalog-entry counts (PAINS A/B/C, BRENK)
  murcko                MurckoScaffoldSmiles (includeChirality=False)

IMPORTANT interpreter note: RDKit 2026.03 changed NumHBA to the strict
acceptor definition (caffeine: 6 → 3). The vendored page wasm is 2026.03.6,
so this fixture MUST be generated with RDKit >= 2026.03 — with 2025.09 the
caffeine HBA reference comes out wrong (6). One known-good interpreter:
  /Library/Frameworks/Python.framework/Versions/3.13/bin/python3 (rdkit 2026.03.6)

The molecule set is chosen for discriminating corners: rhodanine_nmethyl
(PAINS=1), metformin (QED<0.5), cholesterol (Lipinski fail via logP>5,
QED<0.5), hexane (TPSA 0, no Murcko scaffold).

Regenerate with:  python3 scripts/gen_batch_refs.py
"""
import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, FilterCatalog, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

if tuple(int(x) for x in __import__('rdkit').__version__.split('.')[:2]) < (2026, 3):
    sys.exit(
        'RDKit >= 2026.03 required: 2026.03 changed NumHBA to the strict '
        'acceptor definition (caffeine 6 -> 3) and this fixture must match the '
        'vendored 2026.03.6 wasm. Got ' + __import__('rdkit').__version__
    )

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "tests" / "fixtures" / "lbdd" / "batch_refs.json"

MOLECULES = {
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "rhodanine_nmethyl": "CN1C(=O)/C(=C/c2ccccc2)SC1=S",
    "metformin": "CN(C)C(=N)N=C(N)N",
    "cholesterol": "CC(C)CCCC(C)C1CCC2(C)C1CCC1C2CC=C2CC(O)CCC12C",
    "hexane": "CCCCCC",
}


def build_pains_catalog():
    params = FilterCatalog.FilterCatalogParams()
    for c in ("PAINS_A", "PAINS_B", "PAINS_C"):
        params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, c))
    return FilterCatalog.FilterCatalog(params)


def main() -> int:
    pains_cat = build_pains_catalog()
    refs = {}
    bad = []
    for name, smi in MOLECULES.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            bad.append(name)
            continue
        refs[name] = {
            "smiles": Chem.MolToSmiles(m),
            "MW": Descriptors.MolWt(m),
            "cLogP": Crippen.MolLogP(m),
            "TPSA": rdMolDescriptors.CalcTPSA(m),
            "HBD": rdMolDescriptors.CalcNumHBD(m),
            "HBA": rdMolDescriptors.CalcNumHBA(m),
            "RotB": rdMolDescriptors.CalcNumRotatableBonds(m),
            "qed_mean": round(QED.qed(m), 6),
            "qed_max": round(QED.weights_max(m), 6),
            "pains": len(pains_cat.GetMatches(m)),
            "murcko": MurckoScaffold.MurckoScaffoldSmiles(
                smiles=smi, includeChirality=False
            )
            or "",
        }
    if bad:
        print("INVALID SMILES:", bad, file=sys.stderr)
        return 1
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(refs, indent=1))
    print(f"wrote {OUT} ({len(refs)} molecules)")
    for name, v in refs.items():
        print(
            f"  {name:18s} MW {v['MW']:.3f} logP {v['cLogP']:+.3f} TPSA {v['TPSA']:.2f} "
            f"HBD {v['HBD']} HBA {v['HBA']} RotB {v['RotB']} "
            f"QED {v['qed_mean']:.6f} PAINS {v['pains']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
