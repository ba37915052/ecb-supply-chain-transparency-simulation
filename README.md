# ECB Supply Chain Transparency Platform — Simulation Package

Reference implementation of the Monte-Carlo / agent-based computational experiment for:

> Litvin, A. *Conceptual architecture of an industry platform for ensuring the 
> transparency of supply chains for the electronic component base.* [JOURNAL — update on acceptance]

## What this package contains

| File | Contents |
|---|---|
| `simulation.py` | Full experiment: network generation (N=251, M=460, four tiers + OEM), counterfeit introduction/detection model, threshold-cascade with recovery, incident-learning (antifragility), 30-seed replication with 95% CIs and sensitivity analysis |
| `figures.py` | Reproduces all manuscript figures (Fig. X1–X3) |
| `results_main.json` | Summary statistics from 30-replication experiment |
| `results_sensitivity.json` | Parameter sensitivity analysis (p, Δd, θ, shock, ρ) |
| `fig_cpr.png`, `fig_resilience.png`, `fig_sensitivity.png` | Publication-quality figures |
| `requirements.txt` | Python package versions (pinned) |
| `CITATION.cff` | Citation metadata |

## Reproduce the results

```bash
pip install -r requirements.txt
python simulation.py      # ~15–20 min; outputs results_*.json
python figures.py         # creates fig_*.png
```

All randomness is controlled by NumPy PCG64 with fixed seeds, so every run is exactly reproducible while cross-run variability is preserved.

## Headline results

**30 replications × 20,000 shipments each**

| Metric | Baseline (As-Is) | Platform (To-Be) | Effect size |
|---|---|---|---|
| Counterfeit penetration rate (CPR) | 0.0123 | 0.0079 | **−35.7%** [CI 34.5–36.9]; *p* < 0.001; *d* = 8.5 |
| Resilience index (RI) | 0.392 | 0.996 | *p* < 0.001; *d* = −2.19 |
| Cascade size | 0.727 | 0.089 | *p* < 0.001; *d* = 2.52 |
| MTTR (steps) | 36.5 | 1.3 | *p* < 0.001; *d* = 1.70 |

**Sensitivity:** The counterfeit-risk reduction is governed by the platform detection increment Δd (range 18.4%–50.4% for Δd multipliers 0.5–1.5) and is nearly invariant to threat intensity.

## Model assumptions

- Tier cardinalities: 110 (Tier-4) / 75 (Tier-3) / 40 (Tier-2) / 25 (Tier-1) + 1 OEM = 251 nodes
- Recovery probabilities: 0.05 (baseline), 0.20 (platform)
- Learning constant: κ = 0.004
- Cascade threshold: θ = 0.30

All parameters are declared in the manuscript's parameter-provenance table and are covered by sensitivity analysis. See the Limitations section of the manuscript for interpretation constraints.

## License

MIT License. See `LICENSE` for details.

## Citation

If you use this package, please cite both the article and the repository:

```bibtex
@article{Litvin2026,
  author = {Litvin, Alexander},
  title = {Conceptual architecture of an industry platform for ensuring the transparency of supply chains for the electronic component base},
  journal = {[JOURNAL]},
  year = {2026}
}

@software{Litvin2026software,
  author = {Litvin, Alexander},
  title = {ECB Supply Chain Transparency Simulation: Reference Implementation},
  year = {2026},
  url = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

See `CITATION.cff` for machine-readable metadata.
