#!/usr/bin/env python3
"""Cross-model correlation check: WebMM Gaussian shape Tanimoto (fixed pose,
no alignment) vs RDKit hard-sphere grid ShapeTanimotoDist.

These are two different shape models of the same family (analytic Gaussian
overlap vs sphere-occupancy grid), so the expected relation is strong
monotone correlation with a bounded offset — NOT parity. This script
quantifies that relation (Spearman rho, max abs deviation) and prints the
numbers that go into CODE_STATUS.

Usage: python3 scripts/shape_correlation.py
"""
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers


REPO = Path(__file__).resolve().parent.parent
SDFS = sorted((REPO / "tests" / "fixtures" / "conformers").glob("*.sdf"))


def main() -> int:
    pairs = list(itertools.combinations(SDFS, 2))
    ours, grid = [], []
    for a, b in pairs:
        ma = Chem.MolFromMolFile(str(a), removeHs=False, sanitize=True)
        mb = Chem.MolFromMolFile(str(b), removeHs=False, sanitize=True)
        if ma is None or mb is None:
            continue
        grid_sim = 1.0 - rdShapeHelpers.ShapeTanimotoDist(ma, mb)
        out = subprocess.run(
            ["cargo", "run", "--release", "--quiet", "--example", "shape_fixed_pose", "--", str(a), str(b)],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        ours.append(float(out))
        grid.append(grid_sim)
    rho = spearman(ours, grid)
    dev = max(abs(x - y) for x, y in zip(ours, grid))
    print(json.dumps({
        "pairs": len(ours),
        "spearman_rho": round(rho, 4),
        "max_abs_deviation": round(dev, 4),
        "mean_abs_deviation": round(float(np.mean([abs(x - y) for x, y in zip(ours, grid)])), 4),
        "ours_range": [round(min(ours), 4), round(max(ours), 4)],
        "grid_range": [round(min(grid), 4), round(max(grid), 4)],
    }, indent=1))
    return 0


def spearman(x, y):
    """Rank correlation without scipy."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = float(pos)
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


if __name__ == "__main__":
    sys.exit(main())
