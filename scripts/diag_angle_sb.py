"""Extract per-angle and per-SB params from RDKit for glycine and phenol."""

from rdkit import Chem
from rdkit.Chem import AllChem
import math


def get_mol(name, smiles, coords_flat):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    n = len(coords_flat) // 3
    conf = Chem.Conformer(n)
    for i in range(n):
        conf.SetAtomPosition(
            i, (coords_flat[i * 3], coords_flat[i * 3 + 1], coords_flat[i * 3 + 2])
        )
    mol.AddConformer(conf, assignId=True)
    return mol


gly_smiles = "NCC(=O)O"
gly_coords = [
    1.060536,
    -0.359537,
    0.017220,
    1.313786,
    -1.321620,
    -0.316720,
    1.879975,
    0.247839,
    -0.529038,
    -0.418901,
    0.054542,
    0.017610,
    -0.757087,
    -0.435885,
    0.968848,
    -0.780084,
    -0.475278,
    -0.887738,
    -1.607361,
    1.177039,
    0.019348,
    -1.363944,
    2.327917,
    0.036865,
    -2.908414,
    0.805294,
    -0.005805,
    -3.319766,
    1.726115,
    0.013321,
]

phen_smiles = "c1ccc(O)cc1"
phen_coords = [
    -1.583016,
    0.472765,
    0.034105,
    -0.531677,
    1.391140,
    0.036561,
    0.789937,
    0.942578,
    0.006990,
    1.046563,
    -0.425283,
    -0.024844,
    2.321647,
    -0.908022,
    -0.054391,
    0.005052,
    -1.346417,
    -0.027508,
    -1.314635,
    -0.896178,
    0.002067,
    -2.611602,
    0.823582,
    0.057154,
    -0.743749,
    2.456871,
    0.061552,
    1.597221,
    1.667606,
    0.009279,
    2.935484,
    -0.155754,
    -0.048557,
    0.220683,
    -2.410578,
    -0.052520,
    -2.131908,
    -1.612309,
    0.000111,
]


def diag(name, smiles, coords_flat):
    mol = get_mol(name, smiles, coords_flat)
    mp = AllChem.MMFFGetMoleculeProperties(mol, "MMFF94s")
    n = mol.GetNumAtoms()

    print(f"\n{'=' * 70}")
    print(f"  {name}")
    print(f"{'=' * 70}")

    # Atom types
    for i in range(n):
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        at = mp.GetMMFFAtomType(i)
        print(f"  Atom {i}: {sym:2s} -> MMFF type {at}")

    # Find all angles
    conf = mol.GetConformer()
    print(f"\n--- Angle parameters ---")

    total_angle_e = 0.0
    angles = []
    for j in range(n):
        neighbors = [nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors()]
        for idx_a, i in enumerate(neighbors):
            for k in neighbors[idx_a + 1 :]:
                angles.append((i, j, k))

    for i, j, k in angles:
        sym_i = mol.GetAtomWithIdx(i).GetSymbol()
        sym_j = mol.GetAtomWithIdx(j).GetSymbol()
        sym_k = mol.GetAtomWithIdx(k).GetSymbol()
        at_i = mp.GetMMFFAtomType(i)
        at_j = mp.GetMMFFAtomType(j)
        at_k = mp.GetMMFFAtomType(k)

        res = mp.GetMMFFAngleBendParams(mol, i, j, k)
        if res is None:
            print(
                f"  Angle {i}({sym_i})-{j}({sym_j})-{k}({sym_k}): types={at_i}-{at_j}-{at_k} -> NO PARAMS"
            )
            continue

        angle_type, ka, theta0 = res

        # Compute angle
        pi = conf.GetAtomPosition(i)
        pj = conf.GetAtomPosition(j)
        pk = conf.GetAtomPosition(k)
        v1 = [pi.x - pj.x, pi.y - pj.y, pi.z - pj.z]
        v2 = [pk.x - pj.x, pk.y - pj.y, pk.z - pj.z]
        n1 = math.sqrt(sum(x * x for x in v1))
        n2 = math.sqrt(sum(x * x for x in v2))
        dot = sum(a * b for a, b in zip(v1, v2))
        cos_t = max(-1.0, min(1.0, dot / (n1 * n2)))
        theta = math.degrees(math.acos(cos_t))

        # Compute energy: E = 0.5 * c2 * ka * dtheta^2 * (1 + cb * dtheta)
        if theta0 >= 179.0:
            e = 143.9325 * ka * (1.0 + math.cos(math.radians(theta)))
        else:
            dtheta = theta - theta0
            c2 = 143.9325 * (math.pi / 180.0) ** 2
            cb = -0.006981317
            e = 0.5 * c2 * ka * dtheta * dtheta * (1.0 + cb * dtheta)

        total_angle_e += e
        print(
            f"  Angle {i}({sym_i})-{j}({sym_j})-{k}({sym_k}): types={at_i}-{at_j}-{at_k} angleType={angle_type} ka={ka:.6f} theta0={theta0:.3f} theta={theta:.3f} E={e:.6f}"
        )

    print(f"\n  Total angle energy (computed): {total_angle_e:.6f}")

    # Stretch-bend parameters
    print(f"\n--- Stretch-bend parameters ---")
    total_sb_e = 0.0

    for i, j, k in angles:
        sym_i = mol.GetAtomWithIdx(i).GetSymbol()
        sym_j = mol.GetAtomWithIdx(j).GetSymbol()
        sym_k = mol.GetAtomWithIdx(k).GetSymbol()
        at_i = mp.GetMMFFAtomType(i)
        at_j = mp.GetMMFFAtomType(j)
        at_k = mp.GetMMFFAtomType(k)

        res = mp.GetMMFFStretchBendParams(mol, i, j, k)
        if res is None:
            continue

        sb_type, kba_ijk, kba_kji = res

        # Get bond params
        bp_ij = mp.GetMMFFBondStretchParams(mol, i, j)
        bp_kj = mp.GetMMFFBondStretchParams(mol, k, j)
        ap = mp.GetMMFFAngleBendParams(mol, i, j, k)

        if bp_ij and bp_kj and ap:
            _, ka_bond_ij, r0_ij = bp_ij
            _, ka_bond_kj, r0_kj = bp_kj
            _, ka_angle, theta0 = ap

            # Compute bond lengths
            pi = conf.GetAtomPosition(i)
            pj = conf.GetAtomPosition(j)
            pk = conf.GetAtomPosition(k)

            rij = math.sqrt(
                (pi.x - pj.x) ** 2 + (pi.y - pj.y) ** 2 + (pi.z - pj.z) ** 2
            )
            rkj = math.sqrt(
                (pk.x - pj.x) ** 2 + (pk.y - pj.y) ** 2 + (pk.z - pj.z) ** 2
            )

            v1 = [pi.x - pj.x, pi.y - pj.y, pi.z - pj.z]
            v2 = [pk.x - pj.x, pk.y - pj.y, pk.z - pj.z]
            n1 = rij
            n2 = rkj
            dot = sum(a * b for a, b in zip(v1, v2))
            cos_t = max(-1.0, min(1.0, dot / (n1 * n2)))
            theta_rad = math.acos(cos_t)

            dr_ij = rij - r0_ij
            dr_kj = rkj - r0_kj
            dtheta = theta_rad - math.radians(theta0)

            e_sb = 143.9325 * (kba_ijk * dr_ij + kba_kji * dr_kj) * dtheta
            total_sb_e += e_sb

            print(
                f"  SB {i}({sym_i})-{j}({sym_j})-{k}({sym_k}): types={at_i}-{at_j}-{at_k} sbType={sb_type} kbaIJK={kba_ijk:.4f} kbaKJI={kba_kji:.4f} r0ij={r0_ij:.4f} r0kj={r0_kj:.4f} theta0={theta0:.3f} dr_ij={dr_ij:.6f} dr_kj={dr_kj:.6f} dtheta={dtheta:.6f} E={e_sb:.6f}"
            )

    print(f"\n  Total SB energy (computed): {total_sb_e:.6f}")


diag("GLYCINE", gly_smiles, gly_coords)
diag("PHENOL", phen_smiles, phen_coords)
