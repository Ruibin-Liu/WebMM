#!/usr/bin/env python3
"""Audit webmm bond/angle params vs RDKit at full precision across all val sets.

For every bond and angle in every molecule: compare webmm's resolved (kb,r0) /
(ka,theta0) to RDKit's GetMMFFBondStretchParams / GetMMFFAngleBendParams.
Report every divergence > TOL, aggregated by type key.
"""
import json, subprocess, os, sys
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETS = ["val_set", "val_set_new", "val_set_new2", "val_set_new3",
        "val_set_new4", "val_set_new5", "val_set_bulk", "val_set_new6"]
DUMP = os.path.join(ROOT, "target/release/examples/dump_params")
TOL = 1e-4

def rdkit_angles(mol, mp):
    """Enumerate angle triples (i,j_central,k) -> (angle_type, ka, theta0)."""
    nbrs = {a.GetIdx(): [b.GetOtherAtomIdx(a.GetIdx()) for b in a.GetBonds()] for a in mol.GetAtoms()}
    out = {}
    for j in nbrs:
        for ii in range(len(nbrs[j])):
            for kk in range(ii + 1, len(nbrs[j])):
                i, k = nbrs[j][ii], nbrs[j][kk]
                key = (min(i, k), j, max(i, k))
                r = mp.GetMMFFAngleBendParams(mol, i, j, k)
                if r is not None:
                    out[key] = r  # (angle_type, ka, theta0)
    return out

bond_div = {}   # (ti,tj,bo) -> list of (rdkit, webmm)
angle_div = {}  # (ti,tj,tk) -> list of (rdkit, webmm)
nbonds = nangles = 0

for s in SETS:
    d = os.path.join("scripts", s)
    if not os.path.exists(os.path.join(d, "rdkit_ref.json")):
        continue
    web = json.loads(subprocess.run([DUMP, d], capture_output=True, text=True, check=True).stdout)
    for name, w in web.items():
        if "error" in w:
            continue
        sdf = os.path.join(d, f"{name}.sdf")
        mol = Chem.MolFromMolFile(sdf, removeHs=False, sanitize=False)
        if mol is None:
            continue
        mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        if mp is None:
            continue
        # bonds
        wb = {(min(b["i"], b["j"]), max(b["i"], b["j"])): b for b in w["bonds"]}
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            key = (min(i, j), max(i, j))
            if key not in wb:
                continue
            r = mp.GetMMFFBondStretchParams(mol, i, j)
            if r is None:
                continue
            _, kb, r0 = r
            wbnd = wb[key]
            nbonds += 1
            if abs(kb - wbnd["kb"]) > TOL or abs(r0 - wbnd["r0"]) > TOL:
                tk = (wbnd["ti"], wbnd["tj"], wbnd["bo"])
                bond_div.setdefault(tk, []).append(((kb, r0), (wbnd["kb"], wbnd["r0"]), name))
        # angles
        wa = {(min(a["i"], a["k"]), a["j"], max(a["i"], a["k"])): a for a in w["angles"]}
        ra = rdkit_angles(mol, mp)
        for key, (at, ka, t0) in ra.items():
            if key not in wa:
                continue
            nangles += 1
            a = wa[key]
            if abs(ka - a["ka"]) > TOL or abs(t0 - a["t0"]) > TOL:
                tk = (a["ti"], a["tj"], a["tk"])
                angle_div.setdefault(tk, []).append(((ka, t0), (a["ka"], a["t0"]), name))

print(f"Compared {nbonds} bonds, {nangles} angles. Divergent type-keys:")
print(f"\n=== BOND divergences ({len(bond_div)} type-keys) ===")
for tk in sorted(bond_div, key=lambda k: -max(abs(bond_div[k][0][0][0]-bond_div[k][0][1][0]), abs(bond_div[k][0][0][1]-bond_div[k][0][1][1]))):
    r, w, name = bond_div[tk][0]
    print(f"  types{tk[:2]} bo={tk[2]}: rdkit kb={r[0]:.4f} r0={r[1]:.4f} | webmm kb={w[0]:.4f} r0={w[1]:.4f} | dkb={w[0]-r[0]:+.4f} dr0={w[1]-r[1]:+.4f}  ({name}, x{len(bond_div[tk])})")
print(f"\n=== ANGLE divergences ({len(angle_div)} type-keys) ===")
for tk in sorted(angle_div, key=lambda k: -max(abs(angle_div[k][0][0][0]-angle_div[k][0][1][0]), abs(angle_div[k][0][0][1]-angle_div[k][0][1][1])))[:40]:
    r, w, name = angle_div[tk][0]
    print(f"  types{tk}: rdkit ka={r[0]:.4f} t0={r[1]:.4f} | webmm ka={w[0]:.4f} t0={w[1]:.4f} | dka={w[0]-r[0]:+.4f} dt0={w[1]-r[1]:+.4f}  ({name}, x{len(angle_div[tk])})")
print(f"\nTotal: {sum(len(v) for v in bond_div.values())} divergent bonds, {sum(len(v) for v in angle_div.values())} divergent angles")
