#!/usr/bin/env python3
"""
reviewer_runs.py — computations for reviewer comment 6 (sensitivity), comment 9 (control
scenarios) and the Eq. (2) check that supports additional comment a.

The archived model (simulation.py, v1.1.1) is used UNCHANGED: only its existing function
arguments are varied. Common random numbers as in Table 3 (graph 1000+s, shipments 2000+s,
cascades 3000+s, learning 4000+s), so every nominal point reproduces Table 3 exactly; if
results_raw.csv is present, this is asserted automatically.

Usage (next to simulation.py):  python reviewer_runs.py [--n 30] [--out reviewer_results]
Outputs: raw CSVs (one row per replication), summary CSV, SUMMARY.md, fig4_sensitivity.(png|pdf),
analytic_eq2.json, run_manifest.json (versions, hardware, runtime, SHA-256 of code and outputs).
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
import types
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from scipy import stats

import simulation as sim

BASE = dict(theta=0.30, shock_frac=0.10, rho=0.10, p_recover_base=0.05, p_recover_platform=0.20)
SWEEPS = {                                   # comment 6: one-at-a-time around BASE
    "theta": ([0.20, 0.30, 0.40], ("baseline", "platform")),
    "shock_frac": ([0.05, 0.10, 0.20], ("baseline", "platform")),
    "rho": ([0.05, 0.10, 0.15], ("platform",)),               # rho acts in the platform branch only
    "p_recover_platform": ([0.10, 0.15, 0.20, 0.25], ("platform",)),
    "p_recover_base": ([0.05, 0.10], ("baseline",)),
}
CPR_SWEEPS = {"p_scale": [0.5, 1.0, 2.0],            # introduction probabilities p_v
              "d_base_scale": [0.5, 1.0, 1.5],       # baseline detection probabilities d_v
              "delta_scale": [0.5, 1.0, 1.2, 1.5]}   # platform detection increment (protective effect)
ABLATION = {                                 # comment 9: 2x2 protection margin x recovery
    "baseline": dict(platform=False),
    "rho_only": dict(platform=True, rho=0.10, p_recover_platform=0.05),
    "prec_only": dict(platform=True, rho=0.0, p_recover_platform=0.20),
    "both": dict(platform=True, rho=0.10, p_recover_platform=0.20),
}
KAPPAS = [0.0, 0.4]                          # comment 9: learning-rule control
N_SHIP, T_MAX = 20000, 60
TIERS = [4, 3, 2, 1, 0]                      # Tier-4 ... OEM


# ----------------------------------------------------------------------------- statistics
def summary(x, seed=12345):
    x = np.asarray(x, float)
    m, sd, (lo, hi) = sim.ci95(x)
    bs = np.random.default_rng(seed).choice(x, (10000, x.size), replace=True).mean(axis=1)
    return dict(n=x.size, mean=m, sd=sd, ci_low=lo, ci_high=hi,
                boot_low=float(np.percentile(bs, 2.5)), boot_high=float(np.percentile(bs, 97.5)))


def paired(a, b):
    d = np.asarray(a, float) - np.asarray(b, float)
    out = {f"diff_{k}": v for k, v in summary(d).items() if k != "n"}
    out["paired_p"] = float(stats.ttest_1samp(d, 0.0).pvalue) if d.std(ddof=1) > 0 else float("nan")
    return out


# ----------------------------------------------------------------------------- Eq. (2), independent of the simulation code
def eq2_cpr(p, d):
    """Undetected-counterfeit probability for one path: substitution at most once, detection
    only strictly downstream of the substitution node, empty product = 1 (Eq. 2)."""
    total, no_prev = 0.0, 1.0
    for i in range(len(p)):
        surv = np.prod([1 - d[j] for j in range(i + 1, len(d))])
        total += p[i] * no_prev * surv
        no_prev *= 1 - p[i]
    return float(total)


def eq2_mttd(p, d):
    num = den = 0.0
    no_prev = 1.0
    for i in range(len(p)):
        pi, s = p[i] * no_prev, 1.0
        no_prev *= 1 - p[i]
        for j in range(i + 1, len(d)):
            num += pi * s * d[j] * (j - i); den += pi * s * d[j]; s *= 1 - d[j]
    return num / den


def vectors(p_scale=1.0, delta_scale=1.0, platform=False, d_base_scale=1.0):
    p = [sim.P_INTRO[t] * p_scale for t in TIERS[:-1]]
    d = [min(0.99, min(0.99, sim.D_BASE[t] * d_base_scale) + (sim.DELTA_PLATFORM[t] * delta_scale if platform else 0.0))
         for t in TIERS]
    return p, d


def mult(fam, v):
    """(p_scale, delta_scale, d_base_scale) for one point of a CPR sweep."""
    return {"p_scale": (v, 1.0, 1.0), "delta_scale": (1.0, v, 1.0), "d_base_scale": (1.0, 1.0, v)}[fam]


@contextmanager
def scaled_d_base(f):
    """Scale the baseline detection probabilities (module constant D_BASE, read by the model at call
    time) for the duration of a block. The model code is not modified; D_BASE is restored afterwards."""
    orig = dict(sim.D_BASE)
    sim.D_BASE.update({t: min(0.99, v * f) for t, v in orig.items()})
    try:
        yield
    finally:
        sim.D_BASE.clear(); sim.D_BASE.update(orig)


def exact_tie_model():
    """In-memory copy of the archived model in which the failure rule max(0, phi - rho) > theta is
    evaluated as phi > theta + rho (+1e-12), i.e. a node with phi exactly equal to theta + rho does not
    fail, as in exact arithmetic. Used only to quantify the floating-point tie effect
    (0.4 - 0.1 > 0.3 is True in IEEE-754); all reported results use the archived model unchanged."""
    src = Path(sim.__file__).read_text(encoding="utf8")
    old = "if max(0.0, f - eff_rho) > theta:"
    assert src.count(old) == 1, "failure rule not found in simulation.py"
    mod = types.ModuleType("simulation_exact_tie")
    exec(compile(src.replace(old, "if f > theta + eff_rho + 1e-12:"), "simulation_exact_tie", "exec"), mod.__dict__)
    return mod


# ----------------------------------------------------------------------------- runs
def run_cascade(graph, s, platform, **kw):
    G, nodes, oem = graph
    ri, cs, t95 = sim.run_cascade(G, nodes, oem, 3000 + s, platform, **{**BASE, **kw})
    return dict(RI=ri, CS=cs, T95cap=t95, attained=int(t95 < T_MAX))   # attainment at step 60 counted as not attained


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--out", default="reviewer_results")
    a = ap.parse_args()
    n, out = a.n, Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    graphs = [sim.build_graph(1000 + s) for s in range(n)]
    raw = {"cascade_sweeps": [], "ablation_2x2": [], "cpr_sensitivity": [], "learning_control": [], "tie_check": []}

    print("comment 6: cascade sweeps")
    for par, (values, scenarios) in SWEEPS.items():
        for v in values:
            for sc in scenarios:
                for s in range(n):
                    r = run_cascade(graphs[s], s, sc == "platform", **{par: v})
                    raw["cascade_sweeps"].append(dict(family=par, value=v, scenario=sc, rep=s, **r))

    print("comment 9: 2x2 ablation")
    for cell, kw in ABLATION.items():
        kw = dict(kw); is_plat = kw.pop("platform")
        for s in range(n):
            raw["ablation_2x2"].append(dict(cell=cell, rep=s, **run_cascade(graphs[s], s, is_plat, **kw)))

    print("comment 6: CPR sensitivity")
    base_c0 = {}
    for fam, values in CPR_SWEEPS.items():
        for v in values:
            ps, ds, dbs = mult(fam, v)
            with scaled_d_base(dbs):
                for s in range(n):
                    G, nodes, oem = graphs[s]
                    key = (ps, dbs, s)
                    if key not in base_c0:
                        base_c0[key] = sim.run_counterfeit(G, nodes, oem, N_SHIP, 2000 + s, False, ps, 1.0)[0]
                    c1 = sim.run_counterfeit(G, nodes, oem, N_SHIP, 2000 + s, True, ps, ds)[0]
                    c0 = base_c0[key]
                    raw["cpr_sensitivity"].append(dict(family=fam, value=v, rep=s, CPR0=c0, CPR1=c1,
                                                       rel_reduction=(c0 - c1) / c0))

    print("comment 9: learning control")
    for k in KAPPAS:
        for s in range(n):
            G, nodes, oem = graphs[s]
            afi, cprs, _ = sim.run_learning(G, nodes, oem, 4000 + s, kappa=k)
            raw["learning_control"].append(dict(kappa=k, rep=s, AFI=afi, CPR_first=cprs[0],
                                                CPR_last=cprs[-1], CPR_mean=float(np.mean(cprs))))

    print("implementation check: threshold ties")
    exact = exact_tie_model()
    for cell in ("baseline", "rho_only", "both"):
        kw = dict(ABLATION[cell]); is_plat = kw.pop("platform")
        for s in range(n):
            G, nodes, oem = graphs[s]
            ri, cs, t95 = exact.run_cascade(G, nodes, oem, 3000 + s, is_plat, **{**BASE, **kw})
            raw["tie_check"].append(dict(cell=cell, rep=s, RI=ri, CS=cs, T95cap=t95))

    # ------------------------------------------------------------------ summaries
    S = []
    cs_rows = raw["cascade_sweeps"]
    for par, (values, scenarios) in SWEEPS.items():
        for v in values:
            for sc in scenarios:
                sel = [r for r in cs_rows if r["family"] == par and r["value"] == v and r["scenario"] == sc]
                for m in ("RI", "CS", "attained"):
                    S.append(dict(block="cascade_sweep", item=f"{par}={v}", scenario=sc, metric=m,
                                  **summary([r[m] for r in sel])))
            if len(scenarios) == 2:
                b = [r for r in cs_rows if r["family"] == par and r["value"] == v and r["scenario"] == "baseline"]
                p = [r for r in cs_rows if r["family"] == par and r["value"] == v and r["scenario"] == "platform"]
                for m in ("RI", "CS"):
                    S.append(dict(block="cascade_sweep_effect", item=f"{par}={v}", scenario="platform-baseline",
                                  metric=m, **paired([r[m] for r in p], [r[m] for r in b])))
    ab = {c: [r for r in raw["ablation_2x2"] if r["cell"] == c] for c in ABLATION}
    for c in ABLATION:
        for m in ("RI", "CS", "attained"):
            S.append(dict(block="ablation_2x2", item=c, scenario="", metric=m, **summary([r[m] for r in ab[c]])))
    for m in ("RI", "CS"):
        g = {c: np.array([r[m] for r in ab[c]]) for c in ABLATION}
        for name, eff in [("effect_rho", g["rho_only"] - g["baseline"]),
                          ("effect_prec", g["prec_only"] - g["baseline"]),
                          ("effect_both", g["both"] - g["baseline"]),
                          ("interaction", g["both"] - g["rho_only"] - g["prec_only"] + g["baseline"])]:
            S.append(dict(block="ablation_effect", item=name, scenario="", metric=m,
                          **{**summary(eff), "paired_p": float(stats.ttest_1samp(eff, 0).pvalue) if eff.std() > 0 else float("nan")}))
    analytic = {}
    for fam, values in CPR_SWEEPS.items():
        for v in values:
            sel = [r for r in raw["cpr_sensitivity"] if r["family"] == fam and r["value"] == v]
            ps, ds, dbs = mult(fam, v)
            pb, db = vectors(ps, 1.0, False, dbs); pp, dp = vectors(ps, ds, True, dbs)
            an = 1 - eq2_cpr(pp, dp) / eq2_cpr(pb, db)
            sm = summary([r["rel_reduction"] for r in sel])
            S.append(dict(block="cpr_sensitivity", item=f"{fam}={v}", scenario="platform vs baseline",
                          metric="rel_reduction", **sm, analytic_eq2=an,
                          analytic_within_ci=bool(sm["ci_low"] <= an <= sm["ci_high"])))
            analytic[f"{fam}={v}"] = an
    lc = {k: [r for r in raw["learning_control"] if r["kappa"] == k] for k in KAPPAS}
    for k in KAPPAS:
        for m in ("AFI", "CPR_first", "CPR_last"):
            S.append(dict(block="learning_control", item=f"kappa={k}", scenario="", metric=m,
                          **summary([r[m] for r in lc[k]])))
    S.append(dict(block="learning_effect", item="kappa=0.4 minus kappa=0", scenario="", metric="AFI",
                  **paired([r["AFI"] for r in lc[0.4]], [r["AFI"] for r in lc[0.0]])))
    tie = {}
    for cell in ("baseline", "rho_only", "both"):
        arch = [r for r in raw["ablation_2x2"] if r["cell"] == cell]
        ex = [r for r in raw["tie_check"] if r["cell"] == cell]
        tie[cell] = dict(RI_archived=float(np.mean([r["RI"] for r in arch])), RI_exact=float(np.mean([r["RI"] for r in ex])),
                         CS_archived=float(np.mean([r["CS"] for r in arch])), CS_exact=float(np.mean([r["CS"] for r in ex])),
                         reps_changed=sum((a["RI"], a["CS"]) != (e["RI"], e["CS"]) for a, e in zip(arch, ex)))
    b04 = [r["RI"] for r in raw["cascade_sweeps"] if r["family"] == "theta" and r["value"] == 0.40 and r["scenario"] == "baseline"]
    tie["rho_only_exact_equals_baseline_theta_0.40"] = bool(
        all(a == r["RI"] for a, r in zip(b04, [r for r in raw["tie_check"] if r["cell"] == "rho_only"])))

    # ------------------------------------------------------------------ Eq. (2) check (comment a)
    pb, db = vectors(); pp, dp = vectors(platform=True)
    G0, nodes0, oem0 = graphs[0]
    rng = np.random.default_rng(7)
    path_values = []
    for _ in range(5):                       # route independence: Eq. (2) on random feasible paths
        v, path = nodes0[4][rng.integers(len(nodes0[4]))], []
        while True:
            path.append(G0.nodes[v]["tier"])
            if v == oem0:
                break
            succ = list(G0.successors(v)); v = succ[rng.integers(len(succ))]
        p_path = [sim.P_INTRO[t] for t in path[:-1]]
        d_path = [sim.D_BASE[t] for t in path]
        path_values.append(eq2_cpr(p_path, d_path))
    sim_base = [r["CPR0"] for r in raw["cpr_sensitivity"] if r["family"] == "p_scale" and r["value"] == 1.0]
    sim_plat = [r["CPR1"] for r in raw["cpr_sensitivity"] if r["family"] == "p_scale" and r["value"] == 1.0]
    eq2 = dict(CPR_baseline=eq2_cpr(pb, db), CPR_platform=eq2_cpr(pp, dp),
               rel_reduction=1 - eq2_cpr(pp, dp) / eq2_cpr(pb, db),
               MTTD_baseline=eq2_mttd(pb, db), MTTD_platform=eq2_mttd(pp, dp),
               binomial_SD_baseline=float(np.sqrt(eq2_cpr(pb, db) * (1 - eq2_cpr(pb, db)) / N_SHIP)),
               binomial_SD_platform=float(np.sqrt(eq2_cpr(pp, dp) * (1 - eq2_cpr(pp, dp)) / N_SHIP)),
               platform_d_vector_T4_to_OEM=dp,
               route_independence_values=path_values,
               route_independent=bool(np.ptp(path_values) < 1e-15),
               simulated_CPR_baseline=summary(sim_base), simulated_CPR_platform=summary(sim_plat))
    eq2["CPR_baseline_within_ci"] = bool(eq2["simulated_CPR_baseline"]["ci_low"] <= eq2["CPR_baseline"] <= eq2["simulated_CPR_baseline"]["ci_high"])
    eq2["CPR_platform_within_ci"] = bool(eq2["simulated_CPR_platform"]["ci_low"] <= eq2["CPR_platform"] <= eq2["simulated_CPR_platform"]["ci_high"])

    # ------------------------------------------------------------------ self-check against Table 3
    table3 = None
    raw_main = Path(sim.__file__).with_name("results_raw.csv")
    if raw_main.exists() and n == 30:
        main_rows = list(csv.DictReader(raw_main.open()))
        nom_b = [r for r in raw["ablation_2x2"] if r["cell"] == "baseline"]
        nom_p = [r for r in raw["ablation_2x2"] if r["cell"] == "both"]
        exact, ai_rel = [], []
        for s, m in enumerate(main_rows):
            exact += [abs(float(m["RI0"]) - nom_b[s]["RI"]), abs(float(m["RI1"]) - nom_p[s]["RI"]),
                      abs(float(m["CS0"]) - nom_b[s]["CS"]), abs(float(m["CS1"]) - nom_p[s]["CS"]),
                      abs(float(m["CPR0"]) - sim_base[s]), abs(float(m["CPR1"]) - sim_plat[s])]
            ai_rel.append(abs(float(m["AI"]) - lc[0.4][s]["AFI"]) / max(abs(float(m["AI"])), 1e-300))
        # CPR, RI and CS must match exactly on every platform. The AFI slope comes from np.polyfit (LAPACK),
        # whose last-digit rounding differs between BLAS builds (observed on Windows: relative 2.7e-14),
        # so AFI is compared with a relative tolerance of 1e-9.
        table3 = dict(compared_values=len(exact) + len(ai_rel), max_abs_diff_CPR_RI_CS=max(exact),
                      max_rel_diff_AFI=max(ai_rel), reproduces_table3=max(exact) == 0 and max(ai_rel) <= 1e-9)
        assert table3["reproduces_table3"], f"nominal points do not reproduce results_raw.csv: {table3}"

    # ------------------------------------------------------------------ write outputs
    for name, rows in raw.items():
        with (out / f"{name}_raw.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    keys = sorted({k for r in S for k in r}, key=lambda k: ["block", "item", "scenario", "metric"].index(k) if k in ("block", "item", "scenario", "metric") else 9)
    with (out / "summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(S)
    (out / "analytic_eq2.json").write_text(json.dumps(eq2, indent=1))
    (out / "tie_check.json").write_text(json.dumps(tie, indent=1))
    write_markdown(out, S, eq2, table3, n, tie)
    make_figure(out, S, analytic)
    manifest = dict(
        started=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t0)), runtime_s=round(time.time() - t0, 1),
        n_replications=n, seeds=dict(graph="1000+s", shipments="2000+s", cascade="3000+s", learning="4000+s"),
        python=sys.version.split()[0], numpy=np.__version__, scipy=__import__("scipy").__version__,
        networkx=__import__("networkx").__version__, matplotlib=__import__("matplotlib").__version__,
        os=platform.platform(), processor=platform.processor() or platform.machine(), cpu_count=os.cpu_count(),
        ram_gb=ram_gb(), table3_self_check=table3,
        sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [Path(sim.__file__), Path(__file__), *sorted(out.glob("*")) ] if p.is_file() and p.name != "run_manifest.json"})
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"done in {manifest['runtime_s']} s -> {out.resolve()}")


def ram_gb():
    try:
        return round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30, 1)
    except (ValueError, AttributeError, OSError):
        return None


def fmt(x, d=3):
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{d}f}"


def write_markdown(out, S, eq2, table3, n, tie):
    L = [f"# Reviewer comments 6, 9, a — results (n = {n} paired replications)\n"]
    if table3:
        L.append(f"Self-check against results_raw.csv ({table3['compared_values']} values): CPR, RI and CS identical "
                 f"(max |diff| = {table3['max_abs_diff_CPR_RI_CS']}); AFI max relative difference "
                 f"{table3['max_rel_diff_AFI']:.1e} (last-digit rounding of the OLS slope in np.polyfit).\n")
    L += ["## Comment 6 — cascade sensitivity (mean [95% bootstrap CI])", "",
          "| Parameter | Scenario | RI | CS | Share reaching 95% |", "|---|---|---|---|---|"]
    for r in [x for x in S if x["block"] == "cascade_sweep" and x["metric"] == "RI"]:
        cs = next(x for x in S if x["block"] == "cascade_sweep" and x["item"] == r["item"] and x["scenario"] == r["scenario"] and x["metric"] == "CS")
        at = next(x for x in S if x["block"] == "cascade_sweep" and x["item"] == r["item"] and x["scenario"] == r["scenario"] and x["metric"] == "attained")
        L.append(f"| {r['item']} | {r['scenario']} | {fmt(r['mean'])} [{fmt(r['boot_low'])}; {fmt(r['boot_high'])}] | "
                 f"{fmt(cs['mean'])} [{fmt(cs['boot_low'])}; {fmt(cs['boot_high'])}] | {at['mean']:.0%} |")
    L += ["", "## Comment 6 — CPR sensitivity (relative reduction, mean [95% t-CI]; Eq. (2) value)", "",
          "| Multiplier | Simulated | Eq. (2) | Eq. (2) inside CI |", "|---|---|---|---|"]
    for r in [x for x in S if x["block"] == "cpr_sensitivity"]:
        L.append(f"| {r['item']} | {r['mean']:.1%} [{r['ci_low']:.1%}; {r['ci_high']:.1%}] | {r['analytic_eq2']:.1%} | {r['analytic_within_ci']} |")
    L += ["", "## Comment 9 — 2×2 control scenarios (mean [95% bootstrap CI])", "",
          "| Cell | RI | CS | Share reaching 95% |", "|---|---|---|---|"]
    for c in ABLATION:
        g = {m: next(x for x in S if x["block"] == "ablation_2x2" and x["item"] == c and x["metric"] == m) for m in ("RI", "CS", "attained")}
        L.append(f"| {c} | {fmt(g['RI']['mean'])} [{fmt(g['RI']['boot_low'])}; {fmt(g['RI']['boot_high'])}] | "
                 f"{fmt(g['CS']['mean'])} [{fmt(g['CS']['boot_low'])}; {fmt(g['CS']['boot_high'])}] | {g['attained']['mean']:.0%} |")
    L += ["", "| Effect (paired) | ΔRI [95% bootstrap CI] | ΔCS [95% bootstrap CI] |", "|---|---|---|"]
    for e in ("effect_rho", "effect_prec", "effect_both", "interaction"):
        ri = next(x for x in S if x["block"] == "ablation_effect" and x["item"] == e and x["metric"] == "RI")
        cs = next(x for x in S if x["block"] == "ablation_effect" and x["item"] == e and x["metric"] == "CS")
        L.append(f"| {e} | {fmt(ri['mean'])} [{fmt(ri['boot_low'])}; {fmt(ri['boot_high'])}] | {fmt(cs['mean'])} [{fmt(cs['boot_low'])}; {fmt(cs['boot_high'])}] |")
    L += ["", "## Comment 9 — learning-rule control", "", "| κ | AFI mean [95% t-CI] | CPR first batch | CPR last batch |", "|---|---|---|---|"]
    for k in KAPPAS:
        g = {m: next(x for x in S if x["block"] == "learning_control" and x["item"] == f"kappa={k}" and x["metric"] == m) for m in ("AFI", "CPR_first", "CPR_last")}
        L.append(f"| {k} | {g['AFI']['mean']:.2e} [{g['AFI']['ci_low']:.2e}; {g['AFI']['ci_high']:.2e}] | {g['CPR_first']['mean']:.5f} | {g['CPR_last']['mean']:.5f} |")
    e = next(x for x in S if x["block"] == "learning_effect")
    L.append(f"\nAFI difference κ = 0.4 vs κ = 0: {e['diff_mean']:.2e} [{e['diff_ci_low']:.2e}; {e['diff_ci_high']:.2e}], paired p = {e['paired_p']:.3g}")
    L += ["", "## Comment a — Eq. (2) against the simulation", "",
          f"- Eq. (2): CPR {eq2['CPR_baseline']:.6f} (baseline), {eq2['CPR_platform']:.6f} (platform), reduction {eq2['rel_reduction']:.2%}; "
          f"inside simulated 95% CI: {eq2['CPR_baseline_within_ci']}, {eq2['CPR_platform_within_ci']}.",
          f"- MTTD (Eq. 2 logic): {eq2['MTTD_baseline']:.3f} → {eq2['MTTD_platform']:.3f} hops.",
          f"- Platform detection vector, Tier-4 → OEM: {eq2['platform_d_vector_T4_to_OEM']}; Tier-4 detection is never applied.",
          f"- Route independence on 5 random paths: {eq2['route_independent']}.",
          f"- Binomial SD of CPR per replication: {eq2['binomial_SD_baseline']:.2e} / {eq2['binomial_SD_platform']:.2e}."]
    L += ["", "## Implementation note — ties at the failure threshold", "",
          "In double precision 0.4 − 0.1 > 0.3 evaluates as True, so in the platform branch a node whose failed-neighbour "
          "fraction equals exactly θ + ρ fails, whereas in exact arithmetic it does not. The baseline branch (ρ = 0) is unaffected. "
          "All results above use the archived model; the exact-tie variant below only quantifies the effect.", "",
          "| Cell | RI archived | RI exact ties | CS archived | CS exact ties | Replications changed |", "|---|---|---|---|---|---|"]
    for c in ("baseline", "rho_only", "both"):
        t = tie[c]
        L.append(f"| {c} | {t['RI_archived']:.4f} | {t['RI_exact']:.4f} | {t['CS_archived']:.4f} | {t['CS_exact']:.4f} | {t['reps_changed']} |")
    L.append(f"\nWith exact ties, rho_only reproduces baseline θ = 0.40 replication by replication: "
             f"{tie['rho_only_exact_equals_baseline_theta_0.40']} (ρ acts as a shift of the failure threshold).")
    (out / "SUMMARY.md").write_text("\n".join(L), encoding="utf8")


def make_figure(out, S, analytic):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "Times"],
                         "font.size": 8, "pdf.fonttype": 42})
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.3))
    for fam, mk, lab in (("delta_scale", "o", "Δd multiplier"), ("d_base_scale", "^", "d multiplier"),
                         ("p_scale", "s", "p multiplier")):
        rows = [x for x in S if x["block"] == "cpr_sensitivity" and x["item"].startswith(fam)]
        xs = [float(r["item"].split("=")[1]) for r in rows]
        ys = np.array([r["mean"] for r in rows]) * 100
        err = np.array([[r["mean"] - r["ci_low"], r["ci_high"] - r["mean"]] for r in rows]).T * 100
        ax[0].errorbar(xs, ys, yerr=err, marker=mk, capsize=3, lw=1, label=lab)
        ax[0].plot(xs, [analytic[r["item"]] * 100 for r in rows], ls=":", color="k", lw=0.8)
    ax[0].plot([], [], ls=":", color="k", lw=0.8, label="Eq. (2)")
    ax[0].set(xlabel="Parameter multiplier", ylabel="Relative CPR reduction, %", title="(a) Counterfeit penetration")
    ax[0].legend(fontsize=6.5, frameon=False)
    for k, (fam, xl, title) in enumerate((("theta", "Failure threshold θ", "(b) Threshold"),
                                          ("shock_frac", "Initial shock fraction s", "(c) Shock"))):
        for sc, mk in (("baseline", "o--"), ("platform", "s-")):
            rows = [x for x in S if x["block"] == "cascade_sweep" and x["item"].startswith(fam) and x["scenario"] == sc and x["metric"] == "RI"]
            xs = [float(r["item"].split("=")[1]) for r in rows]
            ys = np.array([r["mean"] for r in rows])
            err = np.array([[r["mean"] - r["boot_low"], r["boot_high"] - r["mean"]] for r in rows]).T
            ax[k + 1].errorbar(xs, ys, yerr=err, fmt=mk, capsize=3, lw=1, ms=4, label=sc.capitalize())
        ax[k + 1].set(xlabel=xl, ylabel="Tier-1–OEM connectivity RI", ylim=(-0.03, 1.05), title=title)
        ax[k + 1].legend(fontsize=6.5, frameon=False)
    for a_ in ax:
        a_.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig4_sensitivity.{ext}", dpi=600, bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
