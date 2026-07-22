"""
Reference implementation of the Monte-Carlo / agent-based simulation described in:
Litvin A., "Conceptual architecture of an industry platform for ensuring the
transparency of supply chains for the electronic component base".

Implements:
  1. Four-tier supply graph generation (Tier-4 -> Tier-3 -> Tier-2 -> Tier-1 -> OEM),
     N = 251 nodes, M = 460 edges.
  2. Counterfeit introduction/detection model (per-shipment path sampling).
  3. Threshold cascade model with platform protection margin rho.
  4. Recovery dynamics (MTTR) and resilience index RI (area under functionality curve).
  5. Incident-learning (antifragility) experiment.
  6. Multi-seed replication, summary statistics, Welch t-test, sensitivity analysis.

Python 3.12, numpy 2.4, networkx 3.6, scipy 1.17.
"""

import json
import numpy as np
import networkx as nx
from scipy import stats

# ----------------------------------------------------------------------------
# 1. Network generation
# ----------------------------------------------------------------------------
TIER_SIZES = {4: 110, 3: 75, 2: 40, 1: 25}   # + 1 OEM node = 251 nodes
N_EDGES_TARGET = 460

def build_graph(seed: int) -> nx.DiGraph:
    """Stratified random tier-to-tier DAG with N=251, M=460.

    Construction principles (stated explicitly for reproducibility):
    - Nodes are partitioned into tiers; edges only connect tier t -> tier t-1
      (Tier-1 -> OEM). No intra-tier or tier-skipping edges.
    - Every node has at least one outgoing edge (guaranteed downstream route)
      and every non-Tier-4 node has at least one incoming edge.
    - Remaining edges (up to M=460) are added uniformly at random between
      adjacent tiers, without duplicates.
    """
    rng = np.random.default_rng(seed)
    G = nx.DiGraph()
    nodes = {}
    nid = 0
    for tier, size in TIER_SIZES.items():
        nodes[tier] = list(range(nid, nid + size))
        for v in nodes[tier]:
            G.add_node(v, tier=tier)
        nid += size
    oem = nid
    G.add_node(oem, tier=0)
    nodes[0] = [oem]

    # mandatory connectivity: each node gets >=1 outgoing edge to next tier
    for tier in (4, 3, 2, 1):
        dst_tier = tier - 1
        for v in nodes[tier]:
            u = rng.choice(nodes[dst_tier])
            G.add_edge(v, int(u))
        # each downstream node gets >=1 incoming edge
        for u in nodes[dst_tier]:
            if G.in_degree(u) == 0:
                v = rng.choice(nodes[tier])
                G.add_edge(int(v), u)

    # fill remaining edges uniformly between adjacent tiers
    while G.number_of_edges() < N_EDGES_TARGET:
        tier = rng.choice([4, 3, 2, 1])
        v = int(rng.choice(nodes[tier]))
        u = int(rng.choice(nodes[tier - 1]))
        if not G.has_edge(v, u):
            G.add_edge(v, u)
    return G, nodes, oem

# ----------------------------------------------------------------------------
# 2. Counterfeit penetration experiment
# ----------------------------------------------------------------------------
P_INTRO = {4: 0.003, 3: 0.005, 2: 0.010, 1: 0.012}          # tier-dependent p_v
D_BASE  = {4: 0.15, 3: 0.15, 2: 0.15, 1: 0.35, 0: 0.45}     # baseline d_v
DELTA_PLATFORM = {4: 0.06, 3: 0.06, 2: 0.06, 1: 0.15, 0: 0.15}

def run_counterfeit(G, nodes, oem, n_shipments, seed, platform: bool,
                    p_scale=1.0, delta_scale=1.0):
    rng = np.random.default_rng(seed)
    succ = {v: list(G.successors(v)) for v in G.nodes}
    d = {t: min(0.99, D_BASE[t] + (DELTA_PLATFORM[t] * delta_scale if platform else 0.0))
         for t in D_BASE}
    n_pen = 0          # counterfeit arrived undetected at OEM
    detect_steps = []  # steps from introduction to detection (MTTD)
    n_intro = 0
    t4 = nodes[4]
    for _ in range(n_shipments):
        v = t4[rng.integers(len(t4))]
        counterfeit = False
        step_intro = None
        step = 0
        while True:
            tier = G.nodes[v]['tier']
            # introduction at current node (once per shipment path)
            if not counterfeit and tier != 0:
                if rng.random() < P_INTRO[tier] * p_scale:
                    counterfeit = True
                    step_intro = step
                    n_intro += 1
            # detection at control point of current node (downstream of intro)
            if counterfeit and step_intro is not None and step > step_intro:
                if rng.random() < d[tier]:
                    detect_steps.append(step - step_intro)
                    counterfeit = False   # lot removed
                    break
            if tier == 0:
                if counterfeit:
                    n_pen += 1
                break
            nxt = succ[v]
            v = nxt[rng.integers(len(nxt))]
            step += 1
    cpr = n_pen / n_shipments
    mttd = float(np.mean(detect_steps)) if detect_steps else float('nan')
    return cpr, mttd, n_intro

# ----------------------------------------------------------------------------
# 3. Cascade / resilience experiment
# ----------------------------------------------------------------------------
def run_cascade(G, nodes, oem, seed, platform: bool,
                theta=0.30, shock_frac=0.10, rho=0.10,
                t_max=60, p_recover_base=0.05, p_recover_platform=0.20):
    """Threshold cascade on the undirected dependency graph.

    - Initial shock: shock_frac of Tier-2 and Tier-3 nodes fail.
    - A node fails when the fraction of failed neighbors exceeds theta.
      Platform scenario reduces the *effective* failed-neighbor fraction by
      rho (stress-redistribution damping via multi-sourcing / rerouting).
    - Failed nodes recover independently with per-step probability
      p_recover (higher in the platform scenario: faster rerouting and
      re-qualification supported by verified provenance data).
    - RI = mean over steps of the fraction of Tier-1 nodes with an
      operational path to the OEM (area under the functionality curve).
    - CS = max fraction of failed nodes; MTTR = first step at which Tier-1->OEM
      connectivity returns to >= 95% of its pre-shock value.
    """
    rng = np.random.default_rng(seed)
    U = G.to_undirected()
    nbrs = {v: list(U.neighbors(v)) for v in U.nodes}
    failed = {v: False for v in U.nodes}
    shock_pool = nodes[2] + nodes[3]
    k = int(round(shock_frac * len(shock_pool)))
    for v in rng.choice(shock_pool, size=k, replace=False):
        failed[int(v)] = True
    p_rec = p_recover_platform if platform else p_recover_base
    eff_rho = rho if platform else 0.0

    def t1_connectivity():
        H = G.subgraph([v for v in G.nodes if not failed[v]])
        if oem not in H:
            return 0.0
        ok = 0
        for v in nodes[1]:
            if v in H and nx.has_path(H, v, oem):
                ok += 1
        return ok / len(nodes[1])

    ri_curve = []
    max_failed = sum(failed.values())
    mttr = None
    for t in range(t_max):
        # cascade propagation
        newly = []
        for v in U.nodes:
            if failed[v]:
                continue
            nb = nbrs[v]
            if not nb:
                continue
            f = sum(failed[u] for u in nb) / len(nb)
            if max(0.0, f - eff_rho) > theta:
                newly.append(v)
        for v in newly:
            failed[v] = True
        # peak cascade size is measured AFTER failure propagation and
        # BEFORE recovery (order: fail -> measure peak -> recover -> connectivity)
        max_failed = max(max_failed, sum(failed.values()))
        # recovery
        for v in U.nodes:
            if failed[v] and rng.random() < p_rec:
                failed[v] = False
        c = t1_connectivity()
        ri_curve.append(c)
        if mttr is None and c >= 0.95:
            mttr = t + 1
    ri = float(np.mean(ri_curve))
    cs = max_failed / G.number_of_nodes()
    return ri, cs, (mttr if mttr is not None else t_max)

# ----------------------------------------------------------------------------
# 4. Antifragility (incident learning) experiment
# ----------------------------------------------------------------------------
def run_learning(G, nodes, oem, seed, n_batches=10, batch=5000, kappa=0.4):
    """Platform scenario with incident-learning: after each batch, detection
    probabilities are updated as d <- d + kappa * lambda * (1 - d), where
    lambda = n_detected/batch is the labeled-incident rate of the batch and
    kappa = 0.4 is the learning-rate constant; d is capped at 0.99.
    AI index = -dCPR/dN_incidents (OLS slope of CPR on cumulative incidents)."""
    rng = np.random.default_rng(seed)
    d_extra = {t: DELTA_PLATFORM[t] for t in DELTA_PLATFORM}
    cum_inc, cprs, cums = 0, [], []
    succ = {v: list(G.successors(v)) for v in G.nodes}
    for b in range(n_batches):
        d = {t: min(0.99, D_BASE[t] + d_extra[t]) for t in D_BASE}
        n_pen, n_det = 0, 0
        t4 = nodes[4]
        for _ in range(batch):
            v = t4[rng.integers(len(t4))]
            counterfeit, step_intro, step = False, None, 0
            while True:
                tier = G.nodes[v]['tier']
                if not counterfeit and tier != 0 and rng.random() < P_INTRO[tier]:
                    counterfeit, step_intro = True, step
                if counterfeit and step > step_intro and rng.random() < d[tier]:
                    counterfeit = False
                    n_det += 1
                    break
                if tier == 0:
                    n_pen += counterfeit
                    break
                nxt = succ[v]
                v = nxt[rng.integers(len(nxt))]
                step += 1
        cum_inc += n_det
        cprs.append(n_pen / batch)
        cums.append(cum_inc)
        for t in d_extra:  # learning update
            cur = D_BASE[t] + d_extra[t]
            d_extra[t] = min(0.99 - D_BASE[t], d_extra[t] + kappa * (n_det / batch) * (1 - cur))
    slope = float(np.polyfit(cums, cprs, 1)[0])
    return -slope, cprs, cums

# ----------------------------------------------------------------------------
# 5. Replication experiment
# ----------------------------------------------------------------------------
def ci95(x):
    x = np.asarray(x, float)
    m, s = x.mean(), x.std(ddof=1)
    h = stats.t.ppf(0.975, len(x) - 1) * s / np.sqrt(len(x))
    return m, s, (m - h, m + h)

def main(stage: str = "all"):
    """stage: 'all' (default), 'main' (30-seed experiment only),
    or 'sens' (sensitivity suite only)."""
    N_SEEDS = 30
    N_SHIP = 20000
    if stage in ("all", "main"):
      res = {k: [] for k in ["cpr0", "cpr1", "mttd0", "mttd1",
                              "ri0", "ri1", "cs0", "cs1", "mttr0", "mttr1", "ai"]}
      for s in range(N_SEEDS):
          G, nodes, oem = build_graph(1000 + s)
          c0, m0, _ = run_counterfeit(G, nodes, oem, N_SHIP, 2000 + s, platform=False)
          c1, m1, _ = run_counterfeit(G, nodes, oem, N_SHIP, 2000 + s, platform=True)
          r0, cs0, t0 = run_cascade(G, nodes, oem, 3000 + s, platform=False)
          r1, cs1, t1 = run_cascade(G, nodes, oem, 3000 + s, platform=True)
          ai, _, _ = run_learning(G, nodes, oem, 4000 + s)
          for k, v in zip(res, [c0, c1, m0, m1, r0, r1, cs0, cs1, t0, t1, ai]):
              res[k].append(v)

      out = {}
      for k, v in res.items():
          m, sd, ci = ci95(v)
          out[k] = {"mean": m, "sd": sd, "ci95": ci}
      # relative CPR reduction per seed
      rel = [(a - b) / a for a, b in zip(res["cpr0"], res["cpr1"])]
      m, sd, ci = ci95(rel)
      out["cpr_rel_reduction"] = {"mean": m, "sd": sd, "ci95": ci}
      # Paired t-tests: baseline and platform share the graph and process seeds
      # within each replication, so the design is paired (dz = mean(diff)/sd(diff))
      for a, b, name in [("cpr0", "cpr1", "cpr"), ("mttd0", "mttd1", "mttd"),
                          ("ri0", "ri1", "ri"), ("cs0", "cs1", "cs"),
                          ("mttr0", "mttr1", "mttr")]:
          diff = np.asarray(res[a], float) - np.asarray(res[b], float)
          t, p = stats.ttest_rel(res[a], res[b])
          dz = float(diff.mean() / diff.std(ddof=1))
          out[f"test_{name}"] = {"paired_t": float(t), "p": float(p), "dz": dz}

      # raw per-replication outputs (reviewer requirement 11)
      with open("results_raw.csv", "w") as f:
          f.write("replication,graph_seed,shipment_seed,cascade_seed,learning_seed,"
                  "CPR0,CPR1,MTTD0,MTTD1,RI0,RI1,CS0,CS1,MTTR0,MTTR1,AI\n")
          for s in range(N_SEEDS):
              f.write(f"{s},{1000+s},{2000+s},{3000+s},{4000+s},"
                      f"{res['cpr0'][s]},{res['cpr1'][s]},{res['mttd0'][s]},{res['mttd1'][s]},"
                      f"{res['ri0'][s]},{res['ri1'][s]},{res['cs0'][s]},{res['cs1'][s]},"
                      f"{res['mttr0'][s]},{res['mttr1'][s]},{res['ai'][s]}\n")

      with open("results_main.json", "w") as f:
          json.dump(out, f, indent=1, default=float)
      print(json.dumps(out, indent=1, default=float))

    if stage in ("all", "sens"):
      # ---------------- sensitivity analysis ----------------
      # Each sensitivity replication s uses its own independently generated
      # network (graph seed 1000+s), so intervals reflect BOTH topology and
      # process stochasticity, matching the main-experiment design.
      sens = {}
      base_seeds = range(10)
      graphs = [build_graph(1000 + s) for s in base_seeds]
      def rel_red(p_scale=1.0, delta_scale=1.0):
          r = []
          for s in base_seeds:
              Gs, ns, oe = graphs[s]
              c0, _, _ = run_counterfeit(Gs, ns, oe, N_SHIP, 5000 + s, False, p_scale, 1.0)
              c1, _, _ = run_counterfeit(Gs, ns, oe, N_SHIP, 5000 + s, True, p_scale, delta_scale)
              r.append((c0 - c1) / c0)
          return ci95(r)
      for ps in [0.5, 1.0, 2.0]:
          m, sd, ci = rel_red(p_scale=ps)
          sens[f"p_scale_{ps}"] = {"mean": m, "ci95": ci}
      for ds in [0.5, 1.0, 1.2, 1.5]:
          m, sd, ci = rel_red(delta_scale=ds)
          sens[f"delta_scale_{ds}"] = {"mean": m, "ci95": ci}
      for th in [0.20, 0.30, 0.40]:
          r0s, r1s = [], []
          for s in base_seeds:
              Gs, ns, oe = graphs[s]
              r0, _, _ = run_cascade(Gs, ns, oe, 6000 + s, False, theta=th)
              r1, _, _ = run_cascade(Gs, ns, oe, 6000 + s, True, theta=th)
              r0s.append(r0); r1s.append(r1)
          sens[f"theta_{th}"] = {"ri0": ci95(r0s)[0], "ri1": ci95(r1s)[0]}
      for sf in [0.05, 0.10, 0.20]:
          r0s, r1s = [], []
          for s in base_seeds:
              Gs, ns, oe = graphs[s]
              r0, _, _ = run_cascade(Gs, ns, oe, 7000 + s, False, shock_frac=sf)
              r1, _, _ = run_cascade(Gs, ns, oe, 7000 + s, True, shock_frac=sf)
              r0s.append(r0); r1s.append(r1)
          sens[f"shock_{sf}"] = {"ri0": ci95(r0s)[0], "ri1": ci95(r1s)[0]}
      for rh in [0.05, 0.10, 0.15]:
          r1s = []
          for s in base_seeds:
              Gs, ns, oe = graphs[s]
              r1, _, _ = run_cascade(Gs, ns, oe, 8000 + s, True, rho=rh)
              r1s.append(r1)
          sens[f"rho_{rh}"] = {"ri1": ci95(r1s)[0]}
      # recovery-probability sensitivity (author-defined dynamics parameters)
      for pr in [0.10, 0.15, 0.20, 0.25]:
          r1s, t1s = [], []
          for s in base_seeds:
              Gs, ns, oe = graphs[s]
              r1, _, t1 = run_cascade(Gs, ns, oe, 9000 + s, True,
                                      p_recover_platform=pr)
              r1s.append(r1); t1s.append(t1)
          sens[f"prec_platform_{pr}"] = {"ri1": ci95(r1s)[0],
                                         "ri1_ci95": ci95(r1s)[2],
                                         "mttr1": ci95(t1s)[0]}
      for pr in [0.05, 0.10]:
          r0s = []
          for s in base_seeds:
              Gs, ns, oe = graphs[s]
              r0, _, _ = run_cascade(Gs, ns, oe, 9500 + s, False,
                                     p_recover_base=pr)
              r0s.append(r0)
          sens[f"prec_baseline_{pr}"] = {"ri0": ci95(r0s)[0],
                                         "ri0_ci95": ci95(r0s)[2]}
      with open("results_sensitivity.json", "w") as f:
          json.dump(sens, f, indent=1, default=float)
      print(json.dumps(sens, indent=1, default=float))

if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
