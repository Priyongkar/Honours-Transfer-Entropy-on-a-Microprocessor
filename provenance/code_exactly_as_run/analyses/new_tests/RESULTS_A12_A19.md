# Results of A12–A19

Run 3–4 November 2026 (23:59–01:56). Design: `PREDICTIONS_A12.md`, SHA-256 dc288a0f…9680, logged at the start of
every log. Outputs: `output/A12_results.json` (A12, A13 recording, A14 recording, A18, A19),
`output/A13_A16_chip_results.json`, `output/A15_results.json`, `output/A17_results.json`. All are **hash-logged**; A12, A13
and A14 are partly anticipated (stated in the design). Bootstrap: target clusters, 1,000 draws, seed 2026.
Re-simulation for A16 reproduced the saved chip states with 0 mismatches on both programs.

## Scoreboard

| Test | Prediction | Result | Met? |
|---|---|---|---|
| A12.1 (primary) | AUC(mixed) − AUC(strict) > 0 | +0.155 [0.078, 0.223] | yes |
| A12.2 | strict below chance | 0.449 [0.402, 0.488] | yes |
| A12.3 | mixed above chance | 0.604 [0.536, 0.657] | yes |
| A13.1 (primary) | covered latch pairs below chance, chip | P1 0.228 [0.154, 0.310]; P2 0.103 [0.067, 0.145] | yes, both |
| A13.2 (primary) | uncovered − covered > 0 | P1 +0.200 [0.045, 0.361]; P2 +0.371 [0.267, 0.472]; rec +0.201 [0.076, 0.334] | yes, all three |
| A13.3 | uncovered not below chance | P1 0.428 [0.301, 0.564]; P2 0.474 [0.381, 0.571]; **rec 0.239 [0.118, 0.376]** | chip yes; **recording no** |
| A14.1 (primary) | AUC(B×q₀×π) − AUC(B) > 0, recording | +0.003 [−0.035, 0.041] | **no** |
| A14.2 | B×q₀×π AUC includes 0.5 | 0.454 [0.421, 0.496] | **no** |
| A15.1 (primary) | events, λ=1, below chance | 0.118 [0.109, 0.129] | yes |
| A15.2 (primary) | AUC(λ=1) − AUC(λ=0.5) < 0 | −0.383 [−0.405, −0.360] | yes |
| A15.3 | states − events > 0 at λ=1 | +0.390 [0.368, 0.411] | yes |
| A15.4 | lag-0 CMI given parity above chance, all λ | 0.997–0.999 | yes |
| A16.1 (primary) | > half of latch pairs delayed-only when closed | P1 182/182; P2 213/213 | yes, both |
| A16.2 | lag-1 TE, delayed-only vs none, > 0.5 | P1 0.589 [0.571, 0.606]; P2 0.537 [0.526, 0.547] | yes, both |
| A16.3 | lag-1 − lag-0 CMI > 0 on that class | P1 −0.021 [−0.027, −0.014]; P2 −0.040 [−0.048, −0.032] | **no, both (reversed)** |
| A17.1 (primary) | precision − main run's > 0 | 0.173 vs 0.041; +0.132 [0.074, 0.192] | yes |
| A17.2 | precision − R1's includes 0 | +0.001 [−0.019, 0.022] | yes |
| A18.1 (primary) | block-2 group AUC lower limit ≤ 0.5 | 0.472 [0.452, 0.492] | yes |
| A18.2 | block-2 within-target, strict same-phase, lower ≤ 0.5 | 0.460 [0.390, 0.533] | yes |
| A19.1 (primary) | in-phase share ranks latch sources above others | 0.550 [0.498, 0.597] | **no (narrowly)** |
| A19.2 | median in-phase share: latch > 0.5, strict same-phase < 0.5 | 0.946 vs 0.496 | yes |

Primary predictions: 8 of 10 met (A14.1 and A19.1 not). Secondary: 9 of 11 met (A14.2 and A16.3 not); A13.3 also fails on
the recording.

## Details and what each result does and does not show

**A12 (dose–response).** Targets: 74 mixed, 71 locked, 374 strict. Lag-1 TE: mixed 0.604, locked 0.479 [0.380, 0.584],
strict 0.449. The locked bin is not distinguishable from chance either way, so the gradient is shown between the
end bins, not as a monotone trend through the middle bin. Lag-0 CMI is above chance in every bin: 0.696, 0.799, 0.750.

**A13 (lag 2–7 account).** Covered pairs (k ≥ u − 1) are far below chance on both chip programs and on the recording
(0.038 [0.016, 0.077], 210 pairs); uncovered pairs are higher in all three (55, 115 and 36 pairs). On the chip the
uncovered AUC's interval reaches 0.5 (point estimates 0.43 and 0.47), so "not shown below chance" is the accurate
wording, not "at chance". On the recording the uncovered pairs remain below chance (upper 0.376): history coverage
explains part, not all, of the recording's deficit.

**A14 (between-target remainder).** Not accounted for by the target's quiet-cell terms on the recording (+0.003) or
program 1 (−0.007 [−0.036, 0.022]); partly on program 2 (+0.048 [0.015, 0.080], still 0.428 below chance). A9's exact
within-target account stands; the between-target remainder at matched bound remains unexplained. Failed prediction,
to be reported as such.

**A15 (non-idealised synthetic).** With hidden co-drivers and spurious events, the below-chance ranking is graded in
locking: λ = 0.5, 0.7, 0.85, 1.0 give 0.501, 0.255, 0.119, 0.118. State coding removes the penalty (0.496–0.508 in all
conditions) but does not detect the drivers. Lag-0 CMI given parity is near perfect (0.997–0.999); fast drivers act at
lag 0 by construction, so A15.4 confirms the remedy works where its assumption holds and nothing more.

**A16 (timed clamps).** Every latch pair, on both programs, is delayed-only in all 40 moments of its closed-latch set.
Supplementary check run after the results (`supp_a16_open_closed.py`): every latch pair is immediate in its open-latch
set; base rate of delayed-only among all pairs in a set 3–8%. A1's empty delayed class is therefore explained by
timing: random moments mostly found the latch open. On the full delayed-only class, which by construction now holds
no latch pair (each is immediate in its open set), lag-1 TE is above chance (0.589, 0.537) but lag-0 CMI is higher
again (0.609, 0.577): the A1 secondary prediction fails a second time, reversed, on both programs.

**A17 (IDTxl, clock-aligned).** 237 edges on 100 targets; precision against the netlist 0.173 vs the main run's 0.041
(+0.132 [0.074, 0.192]). 39 of the 41 wired edges are at lag 0 (lag-0 precision 0.222; lag-1 0.033). Precision equals
R1's (+0.001 [−0.019, 0.022]): R1's gain came from admitting lag 0, not from the sample misalignment. Latch precision
rises slightly (0.034 vs 0.027; +0.006 [0.002, 0.012]).

**A18 (source block of two).** Group AUC 0.472, same as lag-1 TE on the same samples (0.469). Within-target on strictly
locked same-phase pairs, block-2 lifts lag-1 TE's 0.400 [0.323, 0.480] to 0.460 [0.390, 0.533]: escaping the bound
removes the significant deficit but does not recover the wiring, as predicted.

**A19 (phase-split local TE).** 40 of 41 latch pairs scored (one had zero lag-1 TE). Latch pairs' TE is almost all
in-phase (median share 0.946) while strictly locked same-phase pairs split evenly (0.496). But within the latch
targets, the in-phase share does not significantly rank latch sources above the other sources (0.550, lower 0.498).

## Summary for the chapters

- Strong confirmations: locking dose–response (A12, A15), latch transfer by intervention (A16.1, exceptionless), history
  coverage on the chip (A13), clock alignment is not R1's mechanism (A17), the block escape without recovery (A18).
- Failed or partial: between-target remainder (A14), lag-1 vs lag-0 on delayed influence (A16.3, second failure),
  recording's uncovered latch pairs still below chance (A13.3), in-phase share ranking (A19.1, narrow).
- Status words: hash-logged throughout; A12–A14 partly anticipated; the A16 open-set check is supplementary.
