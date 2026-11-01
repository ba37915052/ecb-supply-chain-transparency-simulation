# Reproduction check of the archived simulation outputs

Result: **PASS** — 2518 of 2551 values identical; all others within tolerance.

Environment: Python 3.14.4, NumPy 2.4.4, SciPy 1.17.1, NetworkX 3.6.1; Windows-10-10.0.19045-SP0.

| File | Field | Values | Identical | Max relative difference | Tolerance |
|---|---|---|---|---|---|
| results_raw.csv | replication | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | graph_seed | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | shipment_seed | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | cascade_seed | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | learning_seed | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | CPR0 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | CPR1 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | MTTD0 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | MTTD1 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | RI0 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | RI1 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | CS0 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | CS1 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | MTTR0 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | MTTR1 | 30 | 30 | 0.0e+00 | 0e+00 |
| results_raw.csv | AI | 30 | 10 | 2.7e-14 | 1e-09 |
| results_main.json | other statistics | 59 | 49 | 3.5e-15 | 1e-12 |
| results_main.json | AI statistics | 4 | 1 | 1.3e-15 | 1e-09 |
| results_sensitivity_raw.csv | family | 320 | 320 | 0.0e+00 | 0e+00 |
| results_sensitivity_raw.csv | param | 320 | 320 | 0.0e+00 | 0e+00 |
| results_sensitivity_raw.csv | replication | 320 | 320 | 0.0e+00 | 0e+00 |
| results_sensitivity_raw.csv | graph_seed | 320 | 320 | 0.0e+00 | 0e+00 |
| results_sensitivity_raw.csv | metric | 320 | 320 | 0.0e+00 | 0e+00 |
| results_sensitivity_raw.csv | value | 320 | 320 | 0.0e+00 | 0e+00 |
| results_sensitivity.json | other statistics | 88 | 88 | 0.0e+00 | 1e-12 |

simulation.py SHA-256 (LF-normalised): `31a26d78ea0939c6c60f9f5ed555126eff5f94294b84217c76eab7eab1295e2e`