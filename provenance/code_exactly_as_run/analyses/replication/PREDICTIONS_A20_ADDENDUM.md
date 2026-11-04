# A20 addendum: handling clamps that do not settle (4 November 2026, 14:45)

The first run stopped at moment 11 of 80: one clamp left the circuit oscillating, and the simulator raised
"did not converge after 400 iterations". No A20 statistic had been computed; the log shows only progress
lines. Rule added before the rerun, without changing any prediction or class definition in
`PREDICTIONS_A20.md` (SHA-256 d4b14ae0…29e9):

- A clamp that does not settle is recorded as *non-convergent* for that transistor and moment, and counts
  as neither E0 nor E1. The simulator is rebuilt and the moment's snapshot restored before continuing, so
  no state from the failed clamp carries over.
- The number of non-convergent clamps, by transistor and by open or closed moment, is reported.
- A latch transistor's classes use only its convergent moments; one with fewer than 5 convergent closed
  moments is reported but not classified.
