# A20: timed clamps on the recorded program (Donkey Kong), fixed before any clamp is applied

Written on 4 November 2026 (afternoon). The author asked for this experiment to be run in the
background before reading this file; the predictions below were drafted by the AI assistant from the
earlier A16 result and the recording's latch reference, and must be read by the author. The script
prints the SHA-256 of this file, of itself and of the patched simulator source at the start of its log.

## Why

The recording's 41 oriented latch group pairs were defined after the multivariate network's errors had
been seen, and the latch relation was confirmed by intervention only on separate test programs (A16).
The authors' simulator now reproduces the recorded Donkey Kong run exactly, so the recorded pairs can be
perturbed in the recorded program itself.

## Method

- **Simulator.** Sim2600 (Python 3 port in `vendor/Sim2600`), with one addition: a wire can be *forced*
  high or low, so that its group takes the forced value during settling (`mycircuitsimulator.pyx`, the
  `_forced` array). Unforced behaviour is unchanged: the first 2,990 rows still equal the published file in
  every entry. Snapshots copy every state array and the scalar state of the CPU, TIA, PIA and console; a
  restore is checked to reproduce the next three half-clocks exactly before the run.
- **Latch transistors.** As in `recording.py`: a transistor whose gate node gates the clock group (group 6)
  and whose two channel nodes are signal nodes, one *driven* (pull-up or pull-down) and one *storage*
  (neither), with the two nodes' transistors in different groups; it underlies one of the 41 oriented group
  pairs.
- **Moments.** 80 rows drawn uniformly (seed 2026) from the rows that the reduction keeps (row ≡ 0 mod 3)
  between row 4,224 (clock-aligned) and row 401,775 (end of the published recording). At each moment a
  latch transistor is *closed* if it is off and *open* if it is on.
- **Clamp.** Restore the moment; force the driven node to the opposite of its current value; settle;
  record whether the storage node differs from the unclamped state (E₀, same half-clock); hold the force
  for three further TIA half-clocks, one CPU phase, and record whether the storage node differs from the
  unclamped run (E₁).
- **Classes.** A latch transistor is *delayed-only when closed* if E₀ = 0 at every closed moment and E₁ = 1
  at at least one; *immediate when open* if E₀ = 1 at at least one open moment. A recording latch group
  pair is *confirmed* if at least one of its latch transistors is both.

## Predictions

| ID | Prediction | Met if |
|---|---|---|
| **A20.1** | Most recording latch group pairs are delayed-only when closed | More than 20 of the 41 pairs have a transistor delayed-only when closed |
| A20.2 | The same pairs are immediate when open | More than 20 of the 41 have a transistor immediate when open |
| A20.3 | Confirmed pairs (both) | Reported; prediction: more than 20 of 41 |

Reported descriptively: per latch transistor, the share of closed moments with E₀ = 0 and E₁ = 1, and of
open moments with E₀ = 1; latch transistors closed or open at fewer than 5 moments.

Status: hash-logged; exploratory (the reference was defined after the multivariate errors were seen); the
first intervention on the recorded program. A failure would mean that the recording's latch reference
does not name delayed influences in the recorded run, and the recording's latch result would then rest on
the chip programs alone.
