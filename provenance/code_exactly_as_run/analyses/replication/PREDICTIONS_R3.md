# R3: the phase-aware remedy on the two held-out programs, fixed before any statistic is computed

Written on 5 November 2026, before any of the scores below was computed on Pitfall or Space Invaders.
`run_remedy_heldout.py` prints the SHA-256 of this file, of itself and of the shared module at the
start of its log.

## What has and has not been seen

- **Seen.** Everything in `output/R2_*_results.json`: lag-1 TE, lag-0 CMI and layout distance on both
  programs (group AUCs, clock-aligned and whole), the R2 within-target and locking results, and the
  latch results. Those results are reported in the thesis (Table 5.3).
- **Not computed on either program before this file.** Zero-lag mutual information conditioned on the
  sample's phase (`MI|parity`), unconditioned zero-lag MI, the phase-class-only score, and source-group
  size. The remedy score was formed on Donkey Kong after its main results (thesis Section 5.7.2) and
  has never been computed on these programs.
- **Status.** Hash-logged; the *score* is held out, the *programs* are not unseen. Supplementary in the
  thesis's vocabulary: formed after the main results, tested on data it was not formed from.

## Pipeline

Identical to R2 (`run_replication.py`: reduction, clock alignment, groups, wiring reference, layout
distance, target-cluster bootstrap with 1,000 multinomial draws over the units, seed 2026). New scores,
each computed per target over all sources, with no target history (k = 0), as on Donkey Kong
(`analyses/supplementary_recompute/recompute_supplementary_extras.py`):

- `MI|parity`: I(X_t ; Y_t | parity_t), parity_t = t mod 2 counted from the first analysed sample.
- `MI`: I(X_t ; Y_t).
- `phase class only`: 0 for a pair whose units are locked to opposite phases, 1 otherwise (ties).
- `source-group size`: number of transistors in the source group.

Lag-1 TE, lag-0 CMI and layout distance are recomputed by the R2 code, so every paired difference uses
the same bootstrap draws.

## Primary analysis

Clock-aligned samples, group level, unstratified AUC of wired (netlist) against all other ordered group
pairs, 95% target-cluster interval. A prediction is **met** on a program if the interval of the stated
quantity excludes the null in the predicted direction; it is **met ×2** if met on both programs,
**partly** if met on one.

## Predictions (Donkey Kong clock-aligned values in brackets; formed from Donkey Kong)

- **R3.1** `MI|parity` ranks the wiring above chance: lower limit of its group AUC > 0.5. [0.661]
- **R3.2** `MI|parity` minus lag-1 TE > 0. [0.661 against 0.396]
- **R3.3** `MI|parity` minus lag-0 CMI > 0: conditioning on the phase improves on unconditioned
  contemporaneous dependence. [0.661 against 0.561]
- **R3.4** `MI|parity` minus phase class only > 0: the remedy carries wiring beyond the phase relation
  itself. [whole recording: 0.737 against 0.647; not computed clock-aligned]
- **R3.5** Layout distance minus `MI|parity` > 0: the remedy stays below the activity-free baseline.
  [0.928 against 0.661]

## Reported, not predicted

Physical and within-target AUCs of `MI|parity`; group AUC of `MI` and of source-group size; all five
quantities on the whole recording.

## What would count against the account

R3.1 or R3.2 failing on either program would mean the remedy recommended by the thesis's protocol
(check 2, Figure 7.1) does not transfer beyond the program it was formed on. R3.4 failing would mean
the remedy's advantage is the phase relation alone. R3.5 failing would contradict F4.
