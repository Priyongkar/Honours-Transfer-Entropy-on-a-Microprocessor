# Results: A11 (run 3 November 2026, 15:22–15:37)

Design: `PREDICTIONS_A11.md` (approved by the author in advance, before reading it). Hashes of the
design file, `a11_synthetic.py`, `common.py` and `a9_decompose.py` are in `A11.log` and
`run_a11.log`. The only pre-run check simulated one network with a seed outside the test set and
printed timing and an identity check. Raw output: `A11_results.json`.

300 synthetic networks (100 binary event units, 50,000 samples): 5 locking strengths λ × 3 coupling
strengths ρ × 20 networks. Fast (within-sample) edges: 2,400 per condition; delayed (one-sample)
edges: 569–631. Within-target AUCs with 95 % target-cluster intervals (2,000 targets per condition,
1,000 draws).

## Lag-1 TE for fast (within-sample) drivers against same-phase non-drivers

| Locking λ (measured index) | ρ = 0.3 | ρ = 0.6 | ρ = 0.9 |
|---|---|---|---|
| 0.5, none (0.51) | 0.518 [0.499, 0.539] | 0.504 [0.484, 0.524] | 0.467 [0.447, 0.486] |
| 0.7 (0.69–0.70) | 0.410 [0.392, 0.430] | 0.251 [0.238, 0.265] | 0.106 [0.099, 0.116] |
| 0.85 (0.84–0.85) | 0.236 [0.222, 0.252] | 0.114 [0.106, 0.122] | 0.026 [0.023, 0.029] |
| 0.95 (0.95) | 0.208 [0.195, 0.221] | 0.099 [0.091, 0.108] | 0.021 [0.018, 0.025] |
| 1.0, strict (1.00) | 0.184 [0.172, 0.197] | 0.088 [0.081, 0.096] | 0.018 [0.015, 0.021] |
| λ = 1 minus λ = 0.5 | −0.334 [−0.358, −0.311] | −0.416 [−0.436, −0.395] | −0.450 [−0.470, −0.428] |

## Predictions

| Prediction | Result |
|---|---|
| A11.1 (primary): at λ = 1, lag-1 TE ranks fast drivers below non-drivers (ρ = 0.6, 0.9) | **Met**: 0.088 and 0.018 |
| A11.2 (primary): blind without locking (interval includes 0.5), all ρ | **Met for ρ = 0.3 and 0.6; not met for ρ = 0.9**: 0.467 [0.447, 0.486], a small penalty even without locking |
| A11.2 (primary): locking makes it worse (λ = 1 minus λ = 0.5 below 0; ρ = 0.6, 0.9) | **Met**: −0.416 and −0.450 (also −0.334 at ρ = 0.3) |
| A11.3 (remedy): lag-0 CMI given parity above 0.5 in all 15 conditions | **Met**: 0.991–1.000 |
| A11.4 (positive control): lag-1 TE for delayed edges above 0.5 in all 15 conditions | **Met**: 0.997–1.000 |

## Descriptive

- **Diagnostic.** Across the 15 conditions, the mean phase-locking index measured from the events
  alone correlates with the fast-edge lag-1 AUC at Spearman −0.68 (p = 0.005): the more strongly
  units lock, the more misleading lag-1 TE becomes, and locking can be measured without the wiring.
- **Check.** On every λ = 1 network the Remark 3 identity holds for strictly locked same-phase
  pairs; followership κ̃ ranks fast drivers above non-drivers at 0.999–1.000.
- Unconditioned lag-0 CMI ranks fast drivers at 0.990–1.000 in every condition.

## Reading

Outside the 6502, in networks built only to have within-sample coupling and a shared rhythm,
lag-1 transfer entropy is not merely blind to fast coupling: once units lock to the rhythm it ranks
true drivers below non-drivers, and the stronger the locking and the coupling, the further below.
Even moderate locking, 70 % of events in the preferred phase, takes the AUC to 0.25 at ρ = 0.6.
Without locking it is near chance, though at the strongest coupling it is already slightly below
(prediction not met). Contemporaneous dependence conditioned on phase recovers the drivers in every
condition, and lag-1 TE recovers delayed coupling in every condition.

**Limits.** The model is idealised: coupling is clean, there is no noise on the coupling other than
ρ, and the contemporaneous measures are near-perfect, far better than on the recording (lag-0 CMI
group AUC 0.632). A11 shows the direction and the dose–response of the penalty, not its size in
real data. It is synthetic, and its mechanism was derived from the recording.
