# Pass 3: verdict-bearing intervals, recomputed and compared with the thesis

Script: `recompute_supplementary_pass3.py`. Outputs: `output/recomputed_supplementary_pass3.json`,
`output/recompute_pass3.log`, self-test `output/selftest_pass3.json` / `output/recompute_pass3_selftest.log`,
and cached common-history TE matrices in `output/cache_pass3/`. Total runtime about 4 minutes.

## Method

The definitions come from the main notebook (cells 47, 50, 52 and 54) and from methods.tex:
- the structural reference is asserted equal to `analysis_pairwise/structural_counts.npz`;
- δ is read from `thesis_v5/signed_effect.npz`;
- the common-history lag-1 TE (k = 1, 2, 4, 8) and lag-0 CMI (k = 1) were not saved. They were recomputed with `tk.te_matrix`, as in cell 52, at about 30 s each. The estimator reproduces the saved lag-1 TE to within 1e-12.

Four sets of 1,000 bootstrap draws are used:
- **target (main generator):** `tk.bootstrap_targets`, seed 2026, the main run's own generator;
- **target (multinomial):** multinomial, seed 2026, as in pass 1;
- **two-way:** targets multinomial seed 2026, sources an independent multinomial seed 2027;
- **vertex:** one multinomial seed-2026 vector applied in both roles.

The self-test passed. With the main run's generator, every point value and every target-cluster interval that this pass recomputed matches `rq1_discrimination.csv`, `rq1_paired_contrasts.csv`, `rq1_direction.csv` and `rq1_signed_effect.csv` exactly at the printed precision, as do the printed tables.

The two-way and vertex draws from the lost supplementary files cannot be regenerated, because their generator was not recorded. The recomputed two-way and vertex intervals are therefore fresh draws of the same schemes. They are expected to differ from the printed ones in the second or third decimal place.

## 1. tab:m-schemes

| Quantity | Printed | Recomputed | Match |
|---|---|---|---|
| H1 physical AUC, estimate | 0.500 | 0.4996 | yes |
| H1 physical, target | [0.411, 0.596] | [0.411, 0.596] (main gen.); [0.412, 0.597] (multinomial) | yes |
| H1 physical, two-way | [0.399, 0.617] | [0.401, 0.627] | differs: re-drawn, upper +0.010; still includes 0.5 |
| H1 physical, vertex | [0.418, 0.568] | [0.421, 0.572] | differs at 3rd decimal; includes 0.5 |
| H1 group AUC, estimate | 0.469 | 0.4689 | yes |
| H1 group, target | [0.445, 0.493] | [0.445, 0.493] | yes |
| H1 group, two-way | [0.418, 0.528] | [0.423, 0.533] | differs at 3rd decimal; includes 0.5 |
| H1 group, vertex | [0.422, 0.533] | [0.422, 0.532] | yes (within 0.001) |
| H1 Holm p | 1 | 1 (raw 0.881, c = 440) | yes |
| H2 ΔAUC group, estimate | +0.163 | +0.1626 | yes |
| H2 group, target | [+0.142, +0.184] | [+0.142, +0.184] | yes |
| H2 group, two-way | [+0.102, +0.244] | [+0.104, +0.236] | differs (re-drawn); above 0 |
| H2 group, vertex | [+0.099, +0.243] | [+0.105, +0.242] | differs (re-drawn); above 0 |
| H2 ΔAUC physical, estimate | +0.094 | +0.0935 | yes |
| H2 physical, target | [+0.065, +0.128] | [+0.065, +0.128] | yes |
| H2 physical, two-way | [+0.035, +0.214] | [+0.038, +0.216] | differs at 3rd decimal; above 0 |
| H2 physical, vertex | [+0.031, +0.233] | [+0.032, +0.232] | yes (within 0.001) |
| H2 Holm p | 0.012 | 0.012 (raw 2/1001, c = 0) | yes |
| H3 Δacc k=1, est. / target | −0.078 [−0.120, −0.031] | −0.078 [−0.120, −0.031] | yes |
| H3 k=1 two-way / vertex | [−0.174, +0.018] / [−0.173, +0.011] | [−0.183, +0.016] / [−0.183, +0.019] | differs (re-drawn); both include 0 in both versions |
| H3 Δacc k=2, est. / target | −0.141 [−0.198, −0.076] | −0.141 [−0.198, −0.076] | yes |
| H3 k=2 two-way / vertex | [−0.295, −0.002] / [−0.290, −0.008] | [−0.303, −0.005] / [−0.294, −0.006] | differs (re-drawn); both below 0 in both versions |
| H3 Δacc k=4, est. / target | −0.083 [−0.138, −0.026] | −0.083 [−0.138, −0.026] | yes |
| H3 k=4 two-way / vertex | [−0.203, +0.034] / [−0.204, +0.024] | [−0.218, +0.036] / [−0.210, +0.034] | differs (re-drawn); include 0 |
| H3 Δacc k=8, est. / target | −0.003 [−0.050, +0.047] | −0.003 [−0.050, +0.047] | yes |
| H3 k=8 two-way / vertex | [−0.092, +0.093] / [−0.084, +0.087] | [−0.094, +0.098] / [−0.087, +0.094] | differs (re-drawn); include 0 |
| H3 AIS-selected k (not printed in this table) | (tab:direction: +0.026 [−0.015, +0.070]) | +0.026 [−0.015, +0.070]; two-way [−0.047, +0.121]; vertex [−0.045, +0.110] | target yes; others not printed |
| H3 Holm p | 1 | 1 (raw = max over k of the two-sided p = 0.907 at k = 8) | yes |
| H4 Δ% negative group, estimate | +12.9 | +12.94 | yes |
| H4 group, target | [+8.3, +17.2] | [+8.3, +17.2] (main gen.); [+8.0, +17.6] (multinomial) | yes |
| H4 group, two-way | [+1.7, +21.1] | [+2.2, +20.8] | differs (re-drawn); above 0 |
| H4 group, vertex | [+2.6, +21.1] | [+2.6, +20.7] | lower yes, upper −0.4 |
| H4 physical, estimate / target | +3.3 [−9.1, +15.1] | +3.27 [−9.1, +15.1] | yes |
| H4 physical, two-way / vertex | not printed | [−14.7, +24.9] / [−8.2, +14.4] | not printed |
| H4 Holm p | 0.012 | 0.012 (raw 2/1001, c = 0) | yes |
| Holm over H1 to H6 | H1 1, H2 .012, H3 1, H4 .012, H5 .012, H6 .012 | same. H5 and H6 raw p come from the saved fractions of positive draws (pass 1 and R1, both 1.0), not from new draws | yes |

The Holm p follows the caption: 2(c+1)/(B+1), where c counts the draws at or beyond the null on the side opposite the estimate, with target-cluster draws at the deciding level. One aside: the H1 *group* AUC under the target bootstrap has a raw p of 0.014. It is not the deciding level, so it is not in the table.

## 2. tab:app-h3-clock: no clock-group source, 1,172 wired group pairs (all unidirectional)

The accuracy, rate-rule and target-interval columns all match exactly at the printed precision. The group accuracies are 0.538 [0.498, 0.577], 0.461 [0.421, 0.503], 0.531 [0.492, 0.572], 0.629 [0.590, 0.672] and 0.665 [0.627, 0.705]. The physical accuracies are 0.535, 0.490, 0.556, 0.612 and 0.632 with their printed intervals. The rate rule is 0.632 for group and 0.590 for physical. Every target-cluster difference interval in the table is reproduced.

The two-way column was re-drawn:

| Row | Printed two-way | Recomputed two-way | Recomputed vertex |
|---|---|---|---|
| k=1 group | [−0.194, 0.023] | [−0.198, 0.020] | [−0.199, 0.025] |
| k=1 physical | [−0.172, 0.076] | [−0.175, 0.073] | [−0.175, 0.069] |
| k=2 group | [−0.322, −0.001] | [−0.330, −0.004] | [−0.319, −0.003] |
| k=2 physical | [−0.259, 0.076] | [−0.265, 0.083] | [−0.256, 0.071] |
| k=4 group | [−0.232, 0.040] | [−0.237, 0.044] | [−0.229, 0.045] |
| k=4 physical | [−0.178, 0.131] | [−0.183, 0.120] | [−0.181, 0.114] |
| k=8 group | [−0.100, 0.107] | [−0.104, 0.118] | [−0.100, 0.111] |
| k=8 physical | [−0.094, 0.150] | [−0.098, 0.158] | [−0.093, 0.150] |
| AIS group | [−0.043, 0.124] | [−0.054, 0.130] | [−0.051, 0.131] |
| AIS physical | [−0.054, 0.163] | [−0.064, 0.154] | [−0.061, 0.154] |

All recomputed two-way intervals have the same sign pattern as the printed ones. The k = 2 group interval stays just below 0, at an upper limit of −0.004 against the printed −0.001.

## 2b. tab:r1-h4: H4 within parity strata

Every printed value is reproduced exactly:
- the percentages and counts in every cell (for example 86.3 (678), 91.8 (112,609), 81.1 (790), 39.7 (635) and 68.2 (266,830));
- the stratum differences: same parity −5.5 [−9.7, −1.4], opposite −12.5 [−46.4, +34.2], mixed −6.2 [−18.8, +6.8];
- the composition-adjusted difference, −5.6 [−9.8, −1.7] under the target bootstrap;
- the physical overall difference, +3.3 [−9.1, +15.1].

The printed two-way interval of the composition-adjusted difference is [−14.6, +0.9]; the recomputed one is [−14.7, +0.9], and the vertex interval is [−14.7, +0.8]. These match within re-draw noise.

## 3. tab:r1-h2 robustness rows

| Row | Printed | Recomputed | Match |
|---|---|---|---|
| Clock excluded, lag-0 physical AUC | 0.641 [0.598, 0.679] | 0.641 [0.598, 0.679] | yes |
| Clock excluded, lag-0 group AUC | 0.682 [0.656, 0.708] | 0.682 [0.656, 0.708] | yes |
| Clock excluded, within-target | 0.761 | 0.761 | yes |
| Clock excluded, Δ physical | +0.169 [+0.133, +0.206] | +0.169 [+0.133, +0.206] | yes |
| Clock excluded, Δ group | +0.204 [+0.179, +0.231] | +0.204 [+0.179, +0.231] | yes |
| Common k=1, lag-0 physical | 0.585 [0.500, 0.672] | 0.585 [0.500, 0.672] | yes |
| Common k=1, lag-0 group | 0.609 [0.591, 0.627] | 0.609 [0.591, 0.627] | yes |
| Common k=1, within-target | 0.620 | 0.620 | yes |
| Common k=1, Δ physical | +0.082 [+0.057, +0.113] | +0.082 [+0.057, +0.113] | yes |
| Common k=1, Δ group | +0.122 [+0.103, +0.140] | +0.122 [+0.103, +0.140] | yes |
| Node level: 1,645 edges, 1,023 nodes | 1,645 / 1,023 | 1,645 / 1,023 | yes |
| Node AUC lag-1 / lag-0 / MI | 0.469 [0.448, 0.492] / 0.615 [0.594, 0.635] / 0.579 [0.560, 0.598] | 0.469 [0.448, 0.492] / 0.615 [0.594, 0.635] / 0.579 [0.560, 0.598] | yes |
| Node-level Δ | +0.145 [+0.127, +0.165] | +0.145 [+0.127, +0.165] | yes |

Node-level intervals use draws clustered by target node with the main generator, seed 2026.

## Could not recompute / not attempted

- **The two-way and vertex draws as originally made.** Their generator was not recorded. Fresh draws of the same schemes are given instead.
- **The rest of tab:r1-h2.** This covers the AP lifts and the node-level AP lift and within-target AUC. It also covers the rows for parity, parity plus global state, rate strata and clock-aligned data, some of which pass 1 recomputed. None of these was in the scope of this pass. Nothing here suggests they would fail to reproduce.

## Consequence for verdict labels

None changes:
- **H1 robust:** the physical AUC interval includes 0.5 under all three schemes.
- **H2 robust:** the group and physical contrast intervals are above 0 under all three schemes.
- **H3 not supported (robust):** no history has an interval above 0 under any scheme, so the asymmetry never beats the rate rule.
- **H4 levels disagree:** the group interval is above 0 under all schemes, and the physical interval includes 0.
- **Holm p-values:** unchanged.

The re-drawn two-way and vertex limits differ from the printed ones by at most 0.015 in AUC or accuracy units (H3 k=4 two-way lower limit) and by 0.5 points for H4.
