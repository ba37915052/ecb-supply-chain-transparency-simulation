#!/usr/bin/env python3
"""
robustness_tests.py — robustness of the paired baseline-vs-platform comparisons (Table 3) to the
distributional shape of the differences: Wilcoxon signed-rank tests and a paired percentile bootstrap
of the mean difference (100,000 resamples, seed 12345), computed from the archived per-replication
outputs results_raw.csv (30 replications). Differences are baseline minus platform.

Wilcoxon: SciPy defaults — zero differences are discarded; the exact null distribution is used when
there are no zero or tied differences (n <= 50), otherwise the normal approximation without continuity
correction. The method actually used is reported for each metric.

Usage:  python robustness_tests.py [--out robustness_results]
"""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
METRICS = [("CPR", "CPR"), ("MTTD", "MTTD, hops"), ("RI", "RI"), ("CS", "Peak cascade size CS"),
           ("MTTR", "T95cap, steps (column MTTR)")]
B, SEED = 100_000, 12345


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="robustness_results")
    out = HERE / ap.parse_args().out
    out.mkdir(parents=True, exist_ok=True)
    raw = HERE / "results_raw.csv"
    rows = list(csv.DictReader(raw.open()))
    res = {}
    for col, label in METRICS:
        d = np.array([float(r[col + "0"]) - float(r[col + "1"]) for r in rows])
        nz = d[d != 0]
        method = "exact" if nz.size == d.size and np.unique(np.abs(nz)).size == nz.size and d.size <= 50 else "approx"
        w = stats.wilcoxon(d, method=method)
        assert np.isclose(w.pvalue, stats.wilcoxon(d).pvalue, rtol=1e-12, atol=0), "method differs from SciPy default"
        idx = np.random.default_rng(SEED).integers(0, d.size, size=(B, d.size))
        lo, hi = np.percentile(d[idx].mean(axis=1), [2.5, 97.5])
        res[col] = dict(label=label, n=int(d.size), zero_differences=int(d.size - nz.size),
                        positive=int((d > 0).sum()), negative=int((d < 0).sum()), mean_difference=float(d.mean()),
                        boot_ci95=[float(lo), float(hi)], wilcoxon_statistic=float(w.statistic),
                        wilcoxon_p=float(w.pvalue), wilcoxon_method=method)
    meta = dict(bootstrap_resamples=B, seed=SEED, ci="percentile", difference="baseline minus platform",
                python=platform.python_version(), numpy=np.__version__, scipy=__import__("scipy").__version__,
                results_raw_sha256_lf=hashlib.sha256(raw.read_bytes().replace(b"\r\n", b"\n")).hexdigest())
    (out / "robustness.json").write_text(json.dumps(dict(results=res, settings=meta), indent=1), encoding="utf8")
    L = ["# Robustness of the paired comparisons (Table 3)", "",
         f"Paired differences baseline minus platform, n = 30 replications; percentile bootstrap of the mean "
         f"difference, {B:,} resamples, seed {SEED}; Wilcoxon signed-rank test (SciPy {meta['scipy']} defaults).", "",
         "| Metric | Mean difference | 95% bootstrap CI | Wilcoxon p | Method | Zero diffs | +/− |", "|---|---|---|---|---|---|---|"]
    for col, r in res.items():
        L.append(f"| {r['label']} | {r['mean_difference']:.4g} | [{r['boot_ci95'][0]:.4g}; {r['boot_ci95'][1]:.4g}] | "
                 f"{r['wilcoxon_p']:.2e} | {r['wilcoxon_method']} | {r['zero_differences']} | {r['positive']}/{r['negative']} |")
    L.append(f"\nLargest Wilcoxon p: {max(r['wilcoxon_p'] for r in res.values()):.2e}; "
             f"all bootstrap intervals exclude zero: {all(r['boot_ci95'][0] * r['boot_ci95'][1] > 0 for r in res.values())}.")
    (out / "ROBUSTNESS.md").write_text("\n".join(L), encoding="utf8")
    print("\n".join(L[4:]))


if __name__ == "__main__":
    main()
