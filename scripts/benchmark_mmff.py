#!/usr/bin/env python3
"""Benchmark WebMM MMFF94s against the RDKit validation set.

Accuracy (default):
  Regenerates WebMM references with the release `dump_types_energy` example and
  compares them against each set's `rdkit_ref.json`:

    * atom type mismatches          must be 0
    * every |E_webmm - E_rdkit|     must be < 0.01 kcal/mol (documented threshold)
    * every reference molecule      must be present and parse cleanly in WebMM

Speed (default, `--no-speed` to skip):
  Times energy+gradient throughput for WebMM (release `bench_mmff` example) and
  RDKit 2025.09.3 (if importable) on the same SDFs/coords with the same
  protocol (warmup 20 ops, then time until --min-ms of work per molecule).
  Prints pooled us/op and us/atom per set plus the WebMM/RDKit ratio.

Sets covered (230 molecules total):
  val_set 130 | val_set_new 41 | val_set_new2 8 | val_set_new3 6 |
  val_set_new4 6 | val_set_new5 5 | val_set_bulk 32 | val_set_new6 2

Flags:
  --speed-only   skip the accuracy comparison
  --no-speed     skip the speed benchmark (accuracy only)
  --min-ms N     per-molecule timing floor for the speed benchmark (default 30)

Exit code: 0 on full pass, 1 on any accuracy regression.
"""
import argparse
import glob
import json
import math
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETS = ["val_set", "val_set_new", "val_set_new2", "val_set_new3",
        "val_set_new4", "val_set_new5", "val_set_bulk", "val_set_new6"]
TOL = 0.01  # kcal/mol, documented match threshold
WARMUP = 20


def types_of(entry):
    return entry.get("types") or entry.get("atom_types")


def compare_set(dumper, name):
    """Return (report_dict|None, failures_list) for one validation set."""
    ref_path = os.path.join("scripts", name, "rdkit_ref.json")
    if not os.path.exists(ref_path):
        return None, [f"{name}: missing {ref_path}"]
    rd = json.load(open(ref_path))

    proc = subprocess.run([dumper, os.path.join("scripts", name)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        return None, [f"{name}: dump_types_energy exited {proc.returncode}: {proc.stderr[:200]}"]
    us = json.loads(proc.stdout)

    missing = [m for m in rd if m not in us]
    errors = [m for m in us if isinstance(us[m], dict) and "error" in us[m]]
    if missing or errors:
        return None, [f"{name}: missing={missing[:5]} errors={errors[:5]}"]

    atoms = type_mismatch = 0
    charge_l1 = charge_n = 0.0
    energy_pairs = []
    for mol in sorted(rd):
        r, u = rd[mol], us[mol]
        rt, ut = types_of(r), types_of(u)
        if not rt or not ut or len(rt) != len(ut):
            return None, [f"{name}/{mol}: type list mismatch"]
        for a, b in zip(rt, ut):
            atoms += 1
            if int(a) != int(b):
                type_mismatch += 1
        rc, uc = r.get("charges", []), u.get("charges", [])
        for cr, cw in zip(rc, uc):
            charge_l1 += abs(cr - cw)
            charge_n += 1
        if "energy" in r and "energy" in u:
            energy_pairs.append((mol, r["energy"], u["energy"]))

    n = len(energy_pairs)
    over = [p for p in energy_pairs if abs(p[2] - p[1]) > TOL]
    report = {
        "mols": n,
        "atoms": atoms,
        "type_mismatch": type_mismatch,
        "charge_l1": charge_l1 / max(charge_n, 1),
        "r": None, "rmsd": None,
        "over": len(over),
        "worst": sorted(energy_pairs, key=lambda p: -abs(p[2] - p[1]))[:3],
    }
    if n:
        xs = [p[1] for p in energy_pairs]
        ys = [p[2] for p in energy_pairs]
        mx, my = sum(xs) / n, sum(ys) / n
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
        sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n)
        sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n)
        report["r"] = cov / (sx * sy) if sx * sy > 0 else 0.0
        report["rmsd"] = math.sqrt(sum((x - y) ** 2 for x, y in zip(xs, ys)) / n)
    return report, []


def bench_webmm(bench_exe, sdf_dir, min_ms):
    """Run the release bench_mmff example; return {name: stats}."""
    proc = subprocess.run([bench_exe, sdf_dir, "--min-ms", str(min_ms)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"bench_mmff exited {proc.returncode}: {proc.stderr[:300]}")
    return json.loads(proc.stdout)


def bench_rdkit(sdf_dir, min_ms):
    """Time RDKit MMFF94s energy+gradient on the same SDFs (same protocol).
    Returns (stats, None) or (None, error_msg) if RDKit is unavailable/fails."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdForceFieldHelpers
    except ImportError as e:
        return None, f"rdkit import failed: {e}"
    stats = {}
    for path in sorted(glob.glob(os.path.join(sdf_dir, "*.sdf"))):
        name = os.path.basename(path)[:-4]
        try:
            mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
            if mol is None:
                stats[name] = {"error": "parse"}
                continue
            t0 = time.perf_counter()
            props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol,
                                                                  mmffVariant="MMFF94s")
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, props)
            build_us = (time.perf_counter() - t0) * 1e6
            for _ in range(WARMUP):
                ff.CalcEnergy()
                ff.CalcGrad()
            t_start = time.perf_counter()
            ops = 0
            while time.perf_counter() - t_start < min_ms / 1e3:
                ff.CalcEnergy()
                ff.CalcGrad()
                ops += 1
            eg_us = (time.perf_counter() - t_start) * 1e6
            us_per_op = eg_us / max(ops, 1)
            stats[name] = {
                "n_atoms": mol.GetNumAtoms(),
                "build_us": build_us,
                "eg_us": eg_us,
                "eg_ops": ops,
                "us_per_op": us_per_op,
                "us_per_atom": us_per_op / max(mol.GetNumAtoms(), 1),
            }
        except Exception as e:  # noqa: BLE001 - per-molecule, keep going
            stats[name] = {"error": str(e)[:120]}
    return stats, None


def pooled(stats):
    """Pooled means over molecules: sum(eg_us)/sum(eg_ops) and /sum(eg_ops*n_atoms)."""
    eg_us = sum(s.get("eg_us", 0.0) for s in stats.values()
                if isinstance(s, dict) and "eg_us" in s)
    ops = sum(s.get("eg_ops", 0) for s in stats.values()
              if isinstance(s, dict) and "eg_ops" in s)
    atom_ops = sum(s.get("eg_ops", 0) * s.get("n_atoms", 0)
                   for s in stats.values()
                   if isinstance(s, dict) and "eg_ops" in s)
    return eg_us, ops, atom_ops


def run_speed(bench_exe, min_ms):
    print("\nSpeed benchmark: MMFF94s energy+gradient (warmup %d, >= %d ms/mol)" % (WARMUP, min_ms))
    webmm_by_set, rdkit_by_set = {}, {}
    rdkit_err = None
    for s in SETS:
        d = os.path.join("scripts", s)
        webmm_by_set[s] = bench_webmm(bench_exe, d, min_ms)
        st, err = bench_rdkit(d, min_ms)
        if err and rdkit_err is None:
            rdkit_err = err
        rdkit_by_set[s] = st or {}
    if rdkit_err:
        print(f"  (RDKit comparison skipped: {rdkit_err})")
    print(f"{'set':<14} {'atom-ops':>9} {'webmm us/op':>12} {'webmm us/atom':>14} "
          f"{'rdkit us/op':>12} {'rdkit us/atom':>14} {'w/r':>6}")
    tw_us, tw_ops, tw_atom_ops = 0.0, 0, 0
    tr_us, tr_ops, tr_atom_ops = 0.0, 0, 0
    for s in SETS:
        wus, wop, wat = pooled(webmm_by_set[s])
        tw_us += wus; tw_ops += wop; tw_atom_ops += wat
        if rdkit_err is None:
            rus, rop, rat = pooled(rdkit_by_set[s])
            tr_us += rus; tr_ops += rop; tr_atom_ops += rat
        else:
            rus = rop = rat = 0
        w_op = wus / max(wop, 1)
        w_at = wus / max(wat, 1)
        r_op = rus / max(rop, 1) if rus else float("nan")
        r_at = rus / max(rat, 1) if rus else float("nan")
        ratio = w_at / r_at if r_at == r_at and r_at > 0 else float("nan")
        print(f"{s:<14} {wat:>9} {w_op:>12.2f} {w_at:>14.4f} "
              f"{r_op:>12.2f} {r_at:>14.4f} {ratio:>6.2f}")
    w_op = tw_us / max(tw_ops, 1)
    w_at = tw_us / max(tw_atom_ops, 1)
    if rdkit_err is None:
        r_op = tr_us / max(tr_ops, 1)
        r_at = tr_us / max(tr_atom_ops, 1)
        ratio = w_at / r_at
    else:
        r_op = r_at = ratio = float("nan")
    print(f"{'TOTAL':<14} {tw_atom_ops:>9} {w_op:>12.2f} {w_at:>14.4f} "
          f"{r_op:>12.2f} {r_at:>14.4f} {ratio:>6.2f}")
    print("  (pooled: sum(eg_us)/sum(ops) and /sum(ops*atoms); "
          "w/r = WebMM us/atom / RDKit us/atom)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--speed-only", action="store_true",
                        help="skip the accuracy comparison")
    parser.add_argument("--no-speed", action="store_true",
                        help="skip the speed benchmark (accuracy only)")
    parser.add_argument("--min-ms", type=int, default=30,
                        help="per-molecule timing floor in ms (default 30)")
    args = parser.parse_args()

    os.chdir(ROOT)
    subprocess.run(["cargo", "build", "--release", "--example", "dump_types_energy"],
                   check=True, stdout=subprocess.DEVNULL)
    dumper = os.path.join(ROOT, "target", "release", "examples", "dump_types_energy")
    if not args.no_speed:
        subprocess.run(["cargo", "build", "--release", "--example", "bench_mmff"],
                       check=True, stdout=subprocess.DEVNULL)
        bench_exe = os.path.join(ROOT, "target", "release", "examples", "bench_mmff")

    ok = True
    if not args.speed_only:
        failures = []
        reports = {}
        for s in SETS:
            reports[s], fail = compare_set(dumper, s)
            failures += fail

        print(f"{'set':<14} {'mols':>5} {'atoms':>6} {'typemm':>6} {'q|d|':>7} "
              f"{'r':>9} {'rmsd':>8} {'|dE|>0.01':>9}")
        total_m = total_a = total_tm = total_over = 0
        for s in SETS:
            r = reports.get(s)
            if r is None:
                print(f"{s:<14}  FAILED (no data)")
                continue
            total_m += r["mols"]
            total_a += r["atoms"]
            total_tm += r["type_mismatch"]
            total_over += r["over"]
            r_str = f"{r['r']:.6f}" if r["r"] is not None else "  n/a"
            rmsd = f"{r['rmsd']:.4f}" if r["rmsd"] is not None else "  n/a"
            print(f"{s:<14} {r['mols']:>5} {r['atoms']:>6} {r['type_mismatch']:>6} "
                  f"{r['charge_l1']:>7.4f} {r_str:>9} {rmsd:>8} {r['over']:>9}")
            for mol, rr, ww, *_ in r["worst"]:
                print(f"    worst: {mol:<26} rdkit={rr:>10.4f} webmm={ww:>10.4f} "
                      f"d={ww - rr:+.4f}")

        print(f"\nTOTAL: {total_m} molecules / {total_a} atoms, "
              f"{total_tm} type mismatches, {total_over} energy deltas > {TOL}")
        expected_total = sum(
            len(json.load(open(os.path.join("scripts", s, "rdkit_ref.json"))))
            for s in SETS
            if os.path.exists(os.path.join("scripts", s, "rdkit_ref.json"))
        )
        ok = (not failures and total_m == expected_total
              and total_tm == 0 and total_over == 0)
        for f in failures:
            print(f"FAIL: {f}")
        print(f"PASS: {expected_total}/{expected_total} match RDKit <0.01 kcal/mol, 0 type mismatches."
              if ok else "FAIL: regression detected.")

    if not args.no_speed:
        run_speed(bench_exe, args.min_ms)

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
