#!/usr/bin/env python3
"""Generate tests/fixtures/lbdd/refs.json — Python-RDKit gold references for the
Workbench LBDD layer (QED / PAINS / BRENK / Murcko).

Each validation molecule records:
  qed_mean / qed_max / qed_none   QED scores (6 dp)
  props                            the 8 QED descriptor inputs
  pains / brenk                    sorted lists of matched catalog entry ids
  murcko                           MurckoScaffoldSmiles (includeChirality=False)

The JSON is consumed by the Playwright acceptance run (node side) to parity-check
the page's JS implementations. Regenerate with:
  python3 scripts/gen_lbdd_refs.py
"""
import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import FilterCatalog, QED
from rdkit.Chem.Scaffolds import MurckoScaffold

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "tests" / "fixtures" / "lbdd" / "refs.json"

# name -> SMILES. Drug-like regulars + edge cases + positive controls.
VALIDATION_SET = {
    # --- regulars (overlap with repo test molecules / demo presets) ---
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "nicotine": "CN1CCCC1c1ccc(nc1)C",
    "metformin": "CN(C)C(=N)N=C(N)N",
    "salicylic_acid": "OC(=O)c1ccccc1O",
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
    "naphthalene": "c1ccc2ccccc2c1",
    "anthracene": "c1ccc2cc3ccccc3cc2c1",
    "hexane": "CCCCCC",
    "alanine": "CC(N)C(=O)O",
    "threonine": "CC(O)C(N)C(=O)O",
    "glucose_alpha": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "cholesterol": "CC(C)CCCC(C)C1CCC2(C)C1CCC1C2CC=C2CC(O)CCC12C",
    "phenylalanine": "N[C@@H](Cc1ccccc1)C(=O)O",
    "diclofenac": "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    "warfarin": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    "propranolol": "CC(C)NCC(O)COc1cccc2ccccc12",
    "metoprolol": "COCCc1ccc(OCC(O)CNC(C)C)cc1",
    "naproxen": "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
    "theophylline": "Cn1c(=O)c2[nH]cnc2n(C)c1=O",
    "benzocaine": "COC(=O)c1ccc(N)cc1",
    "lidocaine": "CCN(CC)CC(=O)Nc1c(C)cccc1C",
    # --- Murcko edge cases ---
    "ethylbenzene": "CCc1ccccc1",                       # one ring + side chain
    "bibenzyl": "c1ccc(CCc2ccccc2)cc1",                 # two rings + linker
    "spiro_decane": "C1CCC2(CC1)CCCCC2",                # spiro, no aromatics
    "decalin": "C1CCC2CCCCC2C1",                        # fused aliphatic
    "piperazine_linker": "c1ccc(CN2CCNCC2)cc1",         # ring + N linker ring
    "acyclic_only": "CCCCOCCCC",                        # Murcko = empty
    # --- PAINS positive controls ---
    "rhodanine_nmethyl": "CN1C(=O)/C(=C/c2ccccc2)SC1=S",  # pains ene_rhod_A(235)
    "catechol": "c1cc(O)c(O)cc1",                            # pains_b catechol_A(92)
    "hydroquinone": "c1ccc(O)cc1O",
    "coumarin_diOH": "O=c1oc2cc(O)c(O)cc2cc1O-c3ccccc3",     # substituted coumarin
    "aniline_dialkyl": "c1ccc(N(C)C)cc1N(C)C",               # anil di-alk family
    "pyrrole_alkyl": "Cc1cc(C)n(-c2ccccc2)c1C",              # pyrrole alkene family
    "ketoene_A": "CC(=O)C=Cc1ccccc1",                        # ene-one-ish
    "sulfonamide_alkene": "C=CCS(=O)(=O)Nc1ccccc1",
    # --- BRENK / alert controls ---
    "nitro_benzene": "O=[N+]([O-])c1ccccc1",                 # BRENK nitro
    "azo": "c1ccc(N=Nc2ccccc2)cc1",                          # BRENK azo/nitro
    "aldehyde": "CCC(=O)C=O",                                # BRENK aldehyde-ish
    "terminal_alkyne": "CC#CC",                              # BRENK alkyne
    "thiol": "CCS",                                          # BRENK thiol
    "sulfonyl_halide": "CS(=O)(=O)Cl",                       # reactive
    "anhydride": "O=C1OC(=O)C1",                             # BRENK anhydride
    "hydrazine": "CNNC",                                     # BRENK hydrazine
    "cyclohexane": "C1CCCCC1",                               # plain ring, no alerts
    "chloramphenicol": "OCC(NC(=O)C(Cl)Cl)c1ccc([N+](=O)[O-])cc1",
    "isoniazid": "NC(=O)Nc1cccnc1",                          # amide + hydrazide-ish
    "warfarin_enol": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
}


def build_catalogs():
    params = FilterCatalog.FilterCatalogParams()
    for c in ("PAINS_A", "PAINS_B", "PAINS_C", "BRENK"):
        params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, c))
    full = FilterCatalog.FilterCatalog(params)
    p2 = FilterCatalog.FilterCatalogParams()
    for c in ("PAINS_A", "PAINS_B", "PAINS_C"):
        p2.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, c))
    pains_only = FilterCatalog.FilterCatalog(p2)
    return full, pains_only


def main() -> int:
    cat, cat_p = build_catalogs()
    refs = {}
    bad = []
    for name, smi in VALIDATION_SET.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            bad.append(name)
            continue
        props = QED.properties(m)
        pains = sorted(h.GetDescription() for h in cat_p.GetMatches(m))
        brenk = sorted(
            set(h.GetDescription() for h in cat.GetMatches(m)) - set(pains)
        )
        refs[name] = {
            "smiles": smi,
            "qed_mean": round(QED.qed(m), 6),
            "qed_max": round(QED.weights_max(m), 6),
            "qed_none": round(QED.weights_none(m), 6),
            "props": {
                "MW": props.MW,
                "ALOGP": props.ALOGP,
                "HBA": props.HBA,
                "HBD": props.HBD,
                "PSA": props.PSA,
                "ROTB": props.ROTB,
                "AROM": props.AROM,
                "ALERTS": props.ALERTS,
            },
            "pains": pains,
            "brenk": brenk,
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
    n_pains = sum(1 for v in refs.values() if v["pains"])
    n_brenk = sum(1 for v in refs.values() if v["brenk"])
    print(f"wrote {OUT}: {len(refs)} molecules, {n_pains} with PAINS, {n_brenk} with BRENK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
