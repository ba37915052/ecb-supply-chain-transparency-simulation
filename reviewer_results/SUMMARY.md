# Reviewer comments 6, 9, a — results (n = 30 paired replications)

Self-check against results_raw.csv (210 values): CPR, RI and CS identical (max |diff| = 0.0); AFI max relative difference 2.7e-14 (last-digit rounding of the OLS slope in np.polyfit).

## Comment 6 — cascade sensitivity (mean [95% bootstrap CI])

| Parameter | Scenario | RI | CS | Share reaching 95% |
|---|---|---|---|---|
| theta=0.2 | baseline | 0.034 [0.029; 0.039] | 0.999 [0.997; 1.000] | 10% |
| theta=0.2 | platform | 0.628 [0.492; 0.761] | 0.547 [0.407; 0.688] | 73% |
| theta=0.3 | baseline | 0.392 [0.260; 0.535] | 0.743 [0.609; 0.869] | 40% |
| theta=0.3 | platform | 0.996 [0.994; 0.998] | 0.110 [0.102; 0.117] | 100% |
| theta=0.4 | baseline | 0.992 [0.987; 0.996] | 0.108 [0.101; 0.116] | 100% |
| theta=0.4 | platform | 1.000 [1.000; 1.000] | 0.068 [0.065; 0.072] | 100% |
| shock_frac=0.05 | baseline | 0.946 [0.887; 0.981] | 0.148 [0.107; 0.216] | 100% |
| shock_frac=0.05 | platform | 0.999 [0.997; 1.000] | 0.057 [0.052; 0.061] | 100% |
| shock_frac=0.1 | baseline | 0.392 [0.260; 0.535] | 0.743 [0.609; 0.869] | 40% |
| shock_frac=0.1 | platform | 0.996 [0.994; 0.998] | 0.110 [0.102; 0.117] | 100% |
| shock_frac=0.2 | baseline | 0.040 [0.035; 0.044] | 0.997 [0.995; 0.999] | 10% |
| shock_frac=0.2 | platform | 0.973 [0.957; 0.984] | 0.251 [0.237; 0.268] | 97% |
| rho=0.05 | platform | 0.996 [0.994; 0.998] | 0.111 [0.103; 0.119] | 100% |
| rho=0.1 | platform | 0.996 [0.994; 0.998] | 0.110 [0.102; 0.117] | 100% |
| rho=0.15 | platform | 0.999 [0.998; 1.000] | 0.105 [0.098; 0.112] | 100% |
| p_recover_platform=0.1 | platform | 0.989 [0.984; 0.994] | 0.116 [0.108; 0.124] | 100% |
| p_recover_platform=0.15 | platform | 0.995 [0.992; 0.997] | 0.112 [0.104; 0.120] | 100% |
| p_recover_platform=0.2 | platform | 0.996 [0.994; 0.998] | 0.110 [0.102; 0.117] | 100% |
| p_recover_platform=0.25 | platform | 0.998 [0.997; 0.999] | 0.109 [0.102; 0.116] | 100% |
| p_recover_base=0.05 | baseline | 0.392 [0.260; 0.535] | 0.743 [0.609; 0.869] | 40% |
| p_recover_base=0.1 | baseline | 0.424 [0.294; 0.564] | 0.740 [0.603; 0.867] | 60% |

## Comment 6 — CPR sensitivity (relative reduction, mean [95% t-CI]; Eq. (2) value)

| Multiplier | Simulated | Eq. (2) | Eq. (2) inside CI |
|---|---|---|---|
| p_scale=0.5 | 36.1% [34.2%; 38.0%] | 36.2% | True |
| p_scale=1.0 | 35.7% [34.5%; 36.9%] | 36.2% | True |
| p_scale=2.0 | 36.4% [35.5%; 37.4%] | 36.3% | True |
| d_base_scale=0.5 | 29.2% [28.4%; 30.0%] | 28.6% | True |
| d_base_scale=1.0 | 35.7% [34.5%; 36.9%] | 36.2% | True |
| d_base_scale=1.5 | 52.2% [50.3%; 54.1%] | 53.2% | True |
| delta_scale=0.5 | 19.1% [18.1%; 20.2%] | 19.0% | True |
| delta_scale=1.0 | 35.7% [34.5%; 36.9%] | 36.2% | True |
| delta_scale=1.2 | 41.9% [40.7%; 43.1%] | 42.6% | True |
| delta_scale=1.5 | 50.9% [49.8%; 52.0%] | 51.6% | True |

## Comment 9 — 2×2 control scenarios (mean [95% bootstrap CI])

| Cell | RI | CS | Share reaching 95% |
|---|---|---|---|
| baseline | 0.392 [0.260; 0.535] | 0.743 [0.609; 0.869] | 40% |
| rho_only | 0.983 [0.976; 0.990] | 0.120 [0.112; 0.128] | 100% |
| prec_only | 0.628 [0.492; 0.761] | 0.547 [0.407; 0.688] | 73% |
| both | 0.996 [0.994; 0.998] | 0.110 [0.102; 0.117] | 100% |

| Effect (paired) | ΔRI [95% bootstrap CI] | ΔCS [95% bootstrap CI] |
|---|---|---|
| effect_rho | 0.591 [0.448; 0.721] | -0.623 [-0.749; -0.489] |
| effect_prec | 0.236 [0.139; 0.343] | -0.196 [-0.316; -0.088] |
| effect_both | 0.604 [0.461; 0.736] | -0.633 [-0.759; -0.499] |
| interaction | -0.223 [-0.330; -0.126] | 0.186 [0.078; 0.307] |

## Comment 9 — learning-rule control

| κ | AFI mean [95% t-CI] | CPR first batch | CPR last batch |
|---|---|---|---|
| 0.0 | -1.08e-07 [-6.48e-07; 4.33e-07] | 0.00773 | 0.00773 |
| 0.4 | 8.08e-07 [3.57e-07; 1.26e-06] | 0.00773 | 0.00683 |

AFI difference κ = 0.4 vs κ = 0: 9.15e-07 [7.52e-07; 1.08e-06], paired p = 2.98e-12

## Comment a — Eq. (2) against the simulation

- Eq. (2): CPR 0.012318 (baseline), 0.007860 (platform), reduction 36.19%; inside simulated 95% CI: True, True.
- MTTD (Eq. 2 logic): 1.614 → 1.521 hops.
- Platform detection vector, Tier-4 → OEM: [0.21, 0.21, 0.21, 0.5, 0.6]; Tier-4 detection is never applied.
- Route independence on 5 random paths: True.
- Binomial SD of CPR per replication: 7.80e-04 / 6.24e-04.

## Implementation note — ties at the failure threshold

In double precision 0.4 − 0.1 > 0.3 evaluates as True, so in the platform branch a node whose failed-neighbour fraction equals exactly θ + ρ fails, whereas in exact arithmetic it does not. The baseline branch (ρ = 0) is unaffected. All results above use the archived model; the exact-tie variant below only quantifies the effect.

| Cell | RI archived | RI exact ties | CS archived | CS exact ties | Replications changed |
|---|---|---|---|---|---|
| baseline | 0.3919 | 0.3919 | 0.7429 | 0.7429 | 0 |
| rho_only | 0.9829 | 0.9921 | 0.1199 | 0.1085 | 23 |
| both | 0.9962 | 0.9987 | 0.1096 | 0.1049 | 19 |

With exact ties, rho_only reproduces baseline θ = 0.40 replication by replication: True (ρ acts as a shift of the failure threshold).