# R4 results: physical delay against phase relation on the simulated chip

Design: `PREDICTIONS_R4.md` (sha256 39a3c0cc8d3ad515…, printed at the start of the log). Run 5 November
2026; output `output/R4_delay_vs_phase_results.json`, log `output/R4_delay_vs_phase.log`. A first start
(`output/R4_delay_vs_phase_failed_start.log`) stopped before any statistic because the program list
imported an assembler not installed in this environment; the script was changed to read the two saved
program folders directly, as `a13_a16_chip.programs()` lists them.

All four predictions met on both programs (P1; P2):
- R4.1 opposite-phase immediate influence, lag-1 TE: 0.605 [0.594, 0.616]; 0.540 [0.530, 0.552].
- R4.2 lag-1 TE minus lag-0 CMI there: +0.089 [+0.075, +0.104]; +0.072 [+0.058, +0.088].
- R4.3 opposite immediate minus same one-step, lag-1 TE: +0.063 [+0.049, +0.077]; +0.027 [+0.013, +0.042].
- R4.4 mixed pairs, one-step minus immediate, lag-1 TE: +0.075 [+0.055, +0.096]; +0.086 [+0.072, +0.105].

Descriptive: same-phase immediate influence is below chance at lag 1 (0.461; 0.402) and found at lag 0
(0.632; 0.549); same-phase one-step influence is best at lag 2 (0.604; 0.554). On mixed pairs one-step
influence ranks higher than immediate at every lag, so R4.4 does not show lag 1 specifically matched to it.

Reading: lag-1 TE follows the delay in samples, which for locked units the phase relation sets, not the
physical delay. R4.1 and R4.2 were partly anticipated by the recording's opposite-parity wired pairs.
