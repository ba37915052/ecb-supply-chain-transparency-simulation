# ECB Supply Chain Transparency Platform — Simulation Package

Reference implementation of the Monte-Carlo / agent-based computational experiment for:

> Litvin, A. *Conceptual architecture of an industry platform for ensuring the transparency
> of supply chains for the electronic component base.* [JOURNAL — update on acceptance]

**Archived at Zenodo:** concept DOI [10.5281/zenodo.21487060](https://doi.org/10.5281/zenodo.21487060) (always resolves to the latest version).

## What this package contains

| File | Contents |
|---|---|
| `simulation.py` | Full experiment: network generation (N=251, M=460, four tiers + OEM), counterfeit introduction/detection, threshold cascade with recovery (order: fail → measure peak → recover → measure connectivity), incident learning, 30-seed replication with 95% CIs, paired t-tests (t_rel) with d_z, and the sensitivity suite. Run stages: `python simulation.py [all|main|sens]` |
| `figures.py` | Reproduces all simulation figures (Figs. 2–4 of the manuscript) |
| `results_main.json` | Summary statistics of the 30-replication experiment (mean, SD, 95% CI, paired tests) |
| `results_raw.csv` | **Raw per-replication outputs**: one row per replication with all seeds and metric values (CPR, MTTD, RI, CS, MTTR, AI for both scenarios) |
| `results_sensitivity.json` | Sensitivity analysis: p ×0.5–2, Δd ×0.5/1.0/1.2/1.5, θ 0.20–0.40, shock 0.05–0.20, ρ 0.05–0.15, recovery probabilities (platform 0.10–0.25; baseline 0.05/0.10); **each replication uses its own graph seed** |
| `fig_*.png` | Publication figures |
| `requirements.txt`, `CITATION.cff`, `LICENSE` | Pinned versions, citation metadata, MIT license |

## Reproduce the results

```bash
pip install -r requirements.txt
python simulation.py main   # 30-seed experiment → results_main.json, results_raw.csv
python simulation.py sens   # sensitivity suite  → results_sensitivity.json
python figures.py           # → fig_*.png
```

All randomness uses NumPy PCG64 with fixed seeds (graph 1000+s; shipments 2000+s; cascades 3000+s; learning 4000+s; sensitivity 5000–9500+s), so every run is exactly reproducible while cross-run variability is preserved.

## Headline results (v1.1.0; 30 replications × 20,000 shipments)

| Metric | Baseline (As-Is) | Platform (To-Be) | Paired t; d_z |
|---|---|---|---|
| CPR (×10⁻³) | 12.31 ± 0.59 | 7.91 ± 0.44 | t = 46.3; d_z = 8.45 (**−35.7%**, CI 34.5–36.9) |
| Resilience index | 0.39 | 0.996 | t = 8.5; d_z = 1.55 |
| Peak cascade size | 0.74 | 0.11 | t = 9.5; d_z = 1.73 |
| MTTR (steps) | 36.5 | 1.3 | t = 6.6; d_z = 1.21 |

Sensitivity: the CPR reduction is governed by the detection increment Δd — 18.9% (×0.5), 35.5% (×1.0), 41.7% (×1.2), 50.8% (×1.5) — and is insensitive to threat intensity.

## Changelog

- **v1.1.0** — cascade peak measured before recovery (order documented); paired statistics (t_rel, d_z) replacing independent Welch tests; learning update stated as d ← d + κ·λ·(1−d) with κ = 0.4 (numerically identical dynamics, notation corrected); raw per-replication outputs (`results_raw.csv`); sensitivity replications use independent graph seeds; Δd × 1.2 and recovery-probability sweeps added; staged CLI.
- v1.0.x — initial release and extended sensitivity checks.

## License and citation

MIT License. Cite the article and this package (see `CITATION.cff`).
