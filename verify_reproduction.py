#!/usr/bin/env python3
"""
verify_reproduction.py — re-runs the archived simulation.py (stages 'main' and 'sens') in a temporary
directory and compares every regenerated value with the archived outputs:
results_raw.csv, results_main.json, results_sensitivity_raw.csv, results_sensitivity.json.
The archived files are not modified.

Pass criterion
- per-replication outputs (CSV): identical, except the learning slope AI (AFI), which is computed by
  np.polyfit (LAPACK) and is compared with a relative tolerance of 1e-9 (last-digit rounding differs
  between BLAS builds);
- derived statistics (JSON: means, SD, 95% CI, paired t, p, d_z): relative tolerance 1e-12, and 1e-9
  for statistics of AI; the observed maximum differences are reported.

Usage:  python verify_reproduction.py [--out verification_results]
"""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
FILES = {"main": ["results_raw.csv", "results_main.json"],
         "sens": ["results_sensitivity_raw.csv", "results_sensitivity.json"]}


def rel(a, b):
    return 0.0 if a == b else abs(a - b) / max(abs(a), abs(b))


def flatten(obj, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, f"{path}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, f"{path}[{i}]")
    else:
        yield path, obj


def compare_csv(arch, new, tol_cols):
    A, B = list(csv.reader(arch.open())), list(csv.reader(new.open()))
    assert A[0] == B[0] and len(A) == len(B), f"{arch.name}: header or row count differs"
    hdr, res = A[0], {}
    for ra, rb in zip(A[1:], B[1:]):
        for col, a, b in zip(hdr, ra, rb):
            r = res.setdefault(col, dict(n=0, exact=0, max_rel=0.0, tol=tol_cols.get(col, 0.0)))
            r["n"] += 1
            try:
                fa, fb = float(a), float(b)
            except ValueError:                       # text fields (family, param labels)
                assert a == b, f"{arch.name}: {col} differs ({a} vs {b})"
                r["exact"] += 1
                continue
            r["exact"] += fa == fb
            r["max_rel"] = max(r["max_rel"], rel(fa, fb))
    return res


def compare_json(arch, new, tol_default, tol_ai):
    A, B = dict(flatten(json.loads(arch.read_text()))), dict(flatten(json.loads(new.read_text())))
    assert A.keys() == B.keys(), f"{arch.name}: keys differ"
    res = {}
    for k in A:
        grp = "AI statistics" if k.startswith("/ai/") else "other statistics"
        r = res.setdefault(grp, dict(n=0, exact=0, max_rel=0.0, tol=tol_ai if grp == "AI statistics" else tol_default))
        r["n"] += 1; r["exact"] += A[k] == B[k]; r["max_rel"] = max(r["max_rel"], rel(float(A[k]), float(B[k])))
    return res


def sha256_lf(p):
    return hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="verification_results")
    out = HERE / ap.parse_args().out
    out.mkdir(parents=True, exist_ok=True)
    t0, report = time.time(), {}
    with tempfile.TemporaryDirectory() as tmp:
        for stage in ("main", "sens"):
            print(f"re-running simulation.py {stage} ...", flush=True)
            ts = time.time()
            subprocess.run([sys.executable, str(HERE / "simulation.py"), stage], cwd=tmp, check=True,
                           stdout=subprocess.DEVNULL)
            report[f"runtime_{stage}_s"] = round(time.time() - ts, 1)
            for name in FILES[stage]:
                a, b = HERE / name, Path(tmp) / name
                report[name] = (compare_csv(a, b, {"AI": 1e-9}) if name.endswith(".csv")
                                else compare_json(a, b, 1e-12, 1e-9))
    checks = [r for f in FILES.values() for name in f for r in report[name].values()]
    report["pass"] = all(r["max_rel"] <= r["tol"] for r in checks)
    report["values_compared"] = sum(r["n"] for r in checks)
    report["values_identical"] = sum(r["exact"] for r in checks)
    report["environment"] = dict(python=platform.python_version(), os=platform.platform(),
                                 processor=platform.processor() or platform.machine(),
                                 numpy=__import__("numpy").__version__, scipy=__import__("scipy").__version__,
                                 networkx=__import__("networkx").__version__)
    report["sha256_lf_normalised"] = {p: sha256_lf(HERE / p) for p in ["simulation.py", *FILES["main"], *FILES["sens"]]}
    report["runtime_total_s"] = round(time.time() - t0, 1)
    (out / "verification.json").write_text(json.dumps(report, indent=1), encoding="utf8")

    L = ["# Reproduction check of the archived simulation outputs", "",
         f"Result: **{'PASS' if report['pass'] else 'FAIL'}** — {report['values_identical']} of {report['values_compared']} "
         f"values identical; all others within tolerance.", "",
         f"Environment: Python {report['environment']['python']}, NumPy {report['environment']['numpy']}, "
         f"SciPy {report['environment']['scipy']}, NetworkX {report['environment']['networkx']}; {report['environment']['os']}.", "",
         "| File | Field | Values | Identical | Max relative difference | Tolerance |", "|---|---|---|---|---|---|"]
    for f in [*FILES["main"], *FILES["sens"]]:
        for col, r in report[f].items():
            L.append(f"| {f} | {col} | {r['n']} | {r['exact']} | {r['max_rel']:.1e} | {r['tol']:.0e} |")
    L += ["", f"simulation.py SHA-256 (LF-normalised): `{report['sha256_lf_normalised']['simulation.py']}`"]
    (out / "VERIFICATION.md").write_text("\n".join(L), encoding="utf8")
    print(f"{'PASS' if report['pass'] else 'FAIL'}: {report['values_identical']} of {report['values_compared']} values identical, "
          f"all others within tolerance" if report["pass"] else "FAIL: see VERIFICATION.md")
    for f in [*FILES["main"], *FILES["sens"]]:
        worst = max(report[f].items(), key=lambda kv: kv[1]["max_rel"])
        print(f"  {f:28s} max rel. diff {worst[1]['max_rel']:.1e} ({worst[0]})")
    print(f"simulation.py sha256 (LF): {report['sha256_lf_normalised']['simulation.py'][:16]}  | {report['runtime_total_s']} s")
    sys.exit(0 if report["pass"] else 1)


if __name__ == "__main__":
    main()
