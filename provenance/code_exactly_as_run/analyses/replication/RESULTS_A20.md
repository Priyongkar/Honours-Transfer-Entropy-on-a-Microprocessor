# A20 results: timed clamps on the recorded Donkey Kong run (4 November 2026)

Design `PREDICTIONS_A20.md` (SHA-256 d4b14ae0…29e9) and addendum `PREDICTIONS_A20_ADDENDUM.md`
(23f117bd…0d70), both logged at the start of the run. Output `output/A20_results.json`. Status:
hash-logged, exploratory (the reference was defined after the multivariate errors were seen). 80 moments
(rows 5,871 to 397,995), 51 latch transistors underlying the recording's 41 oriented latch group pairs;
each transistor closed at 32–48 moments and open at 27–48. Run time 278 s.

| ID | Prediction | Result | Met |
|---|---|---|---|
| **A20.1** | > 20 of 41 pairs delayed-only when closed | **41 of 41** | yes |
| A20.2 | > 20 of 41 pairs immediate when open | **41 of 41** | yes |
| A20.3 | > 20 of 41 confirmed (both) | **41 of 41** | yes |

Per transistor: all 51 are delayed-only when closed and immediate when open. At closed moments the clamp
changed the storage node one CPU phase later and not at once in a median of 100 % of moments (minimum
97.9 %; no closed moment ever showed an immediate effect). At open moments it changed it at once in 100 %.
Non-convergent clamps: 10 of 4,080, all at open moments, on two transistors (five each); both are still
classified on their convergent moments.

Reading: in the recorded program itself, every one of the recording's 41 latch pairs carries a delayed
influence from its driven to its storage side when the latch is closed, and an immediate one when it is
open. The recording's latch reference therefore names real delayed influences in the recorded run, not
only in separate test programs. As for A16, this confirms the reference, not the measure: that lag-1 TE
detects these pairs is the separate result (0.985).
