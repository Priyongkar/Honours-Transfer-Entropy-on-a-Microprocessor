# Results: A9 and A10 (run 3 November 2026, 14:38–14:54)

Design: `PREDICTIONS_A9.md` (SHA-256 651351957c9b…), fixed and approved before the run.
Scripts: `a9_decompose.py` (1af166fe…), `a10_nmos.py` (e6f2d3be…), `common.py` (df3eb7d8…),
`nmos.py` (0c3b47fb…); hashes are also printed in `A9.log`, `A10.log` and `run_a9_a10.log`.
Pre-run checks used synthetic data and a shapes-only load check; no statistic of either
test was computed before the full run. Not committed to git.
Raw outputs: `A9_results.json`, `A10_results.json`.

## A9 — exact decomposition (hash-logged, exploratory)

**Checks (all pass on all three datasets).** For every strictly locked same-phase pair,
lag-1 TE equals the closed form q₀·F(π, c, a) and the bound equals q₀·G(π, a) to within
6 × 10⁻¹⁶ bits; a = r̃(1 − κ̃)/(1 − c̃) to within 10⁻¹⁴; no monotonicity violation.
For the 36,526 pairs locked at 95 % but not strictly, |TE − formula| has median 0,
95th percentile 0.29 millibits, maximum 39 millibits.

Scope: recording 77,006 pairs, 315 wired (of 770 wired same-phase pairs); chip program 1
410,325 node pairs, 685 wired; chip program 2 473,188 node pairs, 910 wired.

| Test | Recording | Chip program 1 | Chip program 2 | Prediction |
|---|---|---|---|---|
| A9.1 followership κ̃, within target | 0.725 [0.655, 0.795] | 0.686 [0.638, 0.734] | 0.719 [0.672, 0.762] | above 0.5: **met ×3** |
| Activity r̃ alone, within target (desc.) | 0.492 [0.436, 0.555] | 0.547 [0.509, 0.588] | 0.534 [0.498, 0.565] | — |
| Lag-1 TE, within target (desc.; = AUC of a) | 0.400 [0.323, 0.480] | 0.465 [0.423, 0.511] | 0.472 [0.425, 0.510] | — |
| A9.2 TE minus activity, within target | −0.093 [−0.136, −0.055] | −0.082 [−0.106, −0.061] | −0.062 [−0.083, −0.044] | below 0: **met ×3** |
| Lag-1 TE within B deciles (desc.) | 0.452 [0.419, 0.487] | 0.425 [0.400, 0.454] | 0.380 [0.353, 0.409] | — |
| A9.3 target c within B deciles | 0.454 [0.416, 0.494] | 0.478 [0.448, 0.509] | 0.351 [0.321, 0.378] | below 0.5: met, **not met**, met |
| Lag-1 TE within B decile × c quintile | 0.480 [0.446, 0.522] | 0.451 [0.427, 0.475] | 0.446 [0.417, 0.478] | — |
| A9.4 primary: difference from B-decile AUC | +0.028 [−0.014, +0.074] | +0.025 [−0.002, +0.052] | +0.066 [+0.031, +0.099] | above 0: **not met**, **not met**, met |
| A9.4 reading (pre-stated rule) | not accounted for | not accounted for | partly accounted for | |

Descriptive, within B deciles: wired pairs' targets spend more time in the quiet cell
(q₀: 0.549, 0.553, 0.639) and have a smaller in-phase share of it (π: 0.446, 0.437, 0.458).

**Reading.**
- *Within a target the inversion is explained exactly.* Lag-1 TE ranks a target's strictly
  locked same-phase sources by one quantity: the source's change rate in the shared phase
  in samples where the target did not change. Wired sources are as active as unwired ones
  (activity at chance), but their targets follow them (0.69–0.73), which removes exactly
  the samples in which lag-1 TE could see them. Being followed lowers their within-target
  rank by 0.06–0.09, to below chance on the recording (0.400). This holds on the recording
  and on both chip programs.
- *Across targets, at matched bound, the inversion is only partly explained.* The
  pre-registered account (the target's own change propensity c) was not established on the
  recording or program 1 and accounts for part of it on program 2. The residual lies in the
  other target terms (q₀, π), which were reported but not tested.
- The identity is exact for strictly locked pairs and near-exact for 95 %-locked pairs.

## A10 — 50 new NMOS circuits (hash-logged, held out)

Seeds 600–649; 27 circuits have scorable latch input → storage pairs (95 pairs);
481 storage → reader pairs. Circuit-cluster bootstrap, 1,000 draws. 6,000 half-steps
(primary); 20,000 agree.

| Test | 6,000 half-steps | 20,000 | Prediction |
|---|---|---|---|
| A10.1 input → storage, lag-1 TE | 0.779 [0.615, 0.873] | 0.792 [0.617, 0.904] | above 0.5: **met** |
| input → storage, lag-0 CMI (desc.) | 0.850 [0.707, 0.939] | 0.859 [0.728, 0.949] | — |
| A10.2 primary: input → storage, lag 1 minus lag 0 | −0.071 [−0.201, +0.047] | −0.067 [−0.180, +0.035] | above 0: **not met** |
| A10.3 storage → reader, lag 0 minus lag 1 | +0.449 [+0.343, +0.551] | +0.441 [+0.331, +0.538] | above 0: **met** |
| A10.4 storage → reader, lag-1 TE | 0.316 [0.243, 0.405] | 0.328 [0.261, 0.421] | not above 0.5: **met** |

Pooled with A8's ten circuits (descriptive, 6,000): input → storage lag 1 0.797, lag 0
0.839, difference −0.042 [−0.172, +0.069]; reader lag 0 minus lag 1 +0.448.

**Reading.** With 27 informative circuits instead of 7, A8's result stands: in the NMOS
circuits lag-1 TE registers latch transfer but is not specific to lag 1, and the point
estimate now favours lag 0. The post-hoc reading of A8, that the storage → reader relation
is combinational and appears at lag 0, is confirmed out of sample.
