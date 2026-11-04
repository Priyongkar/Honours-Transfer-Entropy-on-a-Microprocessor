# R2 results: replication on Pitfall and Space Invaders (4 November 2026)

Design: `PREDICTIONS_R2.md` (SHA-256 062360e8…5a5b, logged at the start of both runs). Outputs:
`output/R2_Pitfall_results.json`, `output/R2_SpaceInvaders_results.json`. Status: hash-logged, held out.
Pipeline check: on the regenerated Donkey Kong recording it reproduces the main run's 519 groups,
histories and lag-1 TE (max difference 3.2e-16 bits) and the clock-aligned start at sample 1,408.

Recordings: 639,997 rows each after cropping; 213,333 samples per phase (Donkey Kong: 133,926, its
published file being truncated); clock-aligned from sample 1,720 (Pitfall) and 1,408 (Space Invaders);
537 and 540 groups; clock group changes in 99.2 % and 99.3 % of samples; 1,523 and 1,531 wired group pairs.

## Scoreboard (clock-aligned; every prediction met on both games)

| ID | Prediction | Pitfall | Space Invaders | Donkey Kong (main) | Met |
|---|---|---|---|---|---|
| **R2.1** | Lag-1 TE group AUC below chance | 0.397 [0.371, 0.421] | 0.399 [0.376, 0.422] | 0.396 [0.372, 0.421] | yes ×2 |
| R2.1b | Physical AUC not above chance | 0.438 [0.351, 0.530] | 0.435 [0.353, 0.535] | 0.440 [0.352, 0.530] | yes ×2 |
| **R2.2** | Lag-0 CMI minus lag-1 TE > 0 | +0.176 [+0.153, +0.198] | +0.188 [+0.164, +0.211] | +0.165 (aligned) | yes ×2 |
| **R2.3** | Most groups phase-locked | 461 of 537 (85.8 %); 418 strict | 463 of 540 (85.7 %); 409 strict | 445 of 519 (85.7 %) | yes ×2 |
| **R2.4** | Followership of wired sources > 0.5; lag-1 TE minus activity < 0 | 0.705 [0.644, 0.767]; −0.041 [−0.062, −0.025] | 0.709 [0.647, 0.764]; −0.043 [−0.070, −0.025] | 0.725; −0.093 | yes ×2 |
| R2.5 | Mixed minus strict targets > 0 | +0.147 [+0.063, +0.218] | +0.155 [+0.074, +0.226] | +0.155 | yes ×2 |
| **R2.6** | Latch pairs above matched controls | 0.962 [0.889, 0.994] (42 pairs) | 0.967 [0.898, 0.997] (42 pairs) | 0.985 (41) | yes ×2 |
| R2.6b | Lag-0 CMI below lag-1 TE on latch pairs | 0.371 [0.266, 0.476] | 0.329 [0.238, 0.397] | 0.413 | yes ×2 |
| R2.7 | Layout distance above lag-0 CMI | +0.356 [+0.329, +0.386] | +0.343 [+0.316, +0.370] | (0.928 vs 0.632) | yes ×2 |

Whole recording: lag-1 group AUC 0.459 and 0.452; physical 0.483 and 0.473; lag-0 minus lag-1 +0.173 and
+0.185. Proposition 3's identity holds on every strictly locked same-phase pair of both games (max
difference 2.3e-16 bits; 99,112 and 94,036 pairs, 744 and 728 wired).

## What does not replicate in level, and must be reported (none was a written prediction)

1. **Within strictly locked targets, lag-1 TE does not rank wired sources below chance** on the new games:
   0.564 [0.507, 0.607] and 0.581 [0.520, 0.628], against 0.400 on Donkey Kong. The penalty itself
   replicates (A9.2: wired sources fall below where their activity places them, −0.041 and −0.043, about
   half the Donkey Kong size); the inversion to below chance within a target is specific to Donkey Kong.
   The group-level inversion (R2.1) replicates in full.
2. **The strict bin of the dose–response is not below chance** on either game (upper limits 0.513 and
   0.512), unlike Donkey Kong (0.449 [0.402, 0.488]); the gradient between the end bins replicates.
3. **Across programs there is no variation to test the descriptive prediction**: locking (85.7–85.8 %),
   lag-1 group AUC (0.396–0.399) and the share of strictly locked same-phase wired pairs (0.48–0.49 on the
   new games) barely differ. The three programs do not separate the account's cross-program prediction.

## Reading

The failure, its group-level size, the phase structure, the followership penalty, the dose–response
between end bins, the latch result, the superiority of lag 0 and of layout distance all replicate on two
programs the account was not derived from, with every written prediction met. The within-target
inversion to below chance does not: it is a property of Donkey Kong, and the general claim is the
penalty relative to activity. Same chip, netlist and simulator: this tests dependence on the program,
not on the system.
