# Robustness of the paired comparisons (Table 3)

Paired differences baseline minus platform, n = 30 replications; percentile bootstrap of the mean difference, 100,000 resamples, seed 12345; Wilcoxon signed-rank test (SciPy 1.17.1 defaults).

| Metric | Mean difference | 95% bootstrap CI | Wilcoxon p | Method | Zero diffs | +/− |
|---|---|---|---|---|---|---|
| CPR | 0.004408 | [0.004233; 0.0046] | 1.73e-06 | approx | 0 | 30/0 |
| MTTD, hops | 0.0954 | [0.0865; 0.1045] | 1.86e-09 | exact | 0 | 30/0 |
| RI | -0.6043 | [-0.738; -0.4635] | 2.56e-06 | approx | 1 | 0/29 |
| Peak cascade size CS | 0.6333 | [0.5016; 0.7595] | 1.72e-06 | approx | 0 | 30/0 |
| T95cap, steps (column MTTR) | 35.2 | [25.23; 45.13] | 4.90e-05 | approx | 11 | 19/0 |

Largest Wilcoxon p: 4.90e-05; all bootstrap intervals exclude zero: True.