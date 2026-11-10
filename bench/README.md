# Benchmarks

`bench/bench.sh [binary] [VOL]` re-finds the four-cover point of the generator of y² = x³ + 730033053750
(`results/k730033053750/cover-c1-156.coef`, u = (2045913748, −1557028403, −858328586, 334741533)) with H = 3·10⁹.

Historical search-kernel timings: AMD EPYC VM (x86_64), gcc 13.3 `-O2`, one thread per process; descent, model reduction and certification are excluded. The 9 min 19 s + 11 min 37 s entry adds two completed half-search process durations. The 12 h entry is approximate time to the point in an unfinished half-search while another half-search occupied a second core. No timing is a certificate of complete coverage.

| program | case | H | VOL | wall time | result |
|---|---|---|---|---|---|
| elkies4q | k = 730033053750, cover c1-156 | 3·10⁹ | 8 | 23 min 48 s | point found |
| elkies4q | same, run as two halves (`F0 F1` = 0 0.5, 0.5 1) | 3·10⁹ | 64 | 9 min 19 s + 11 min 37 s | point found (second half) |
| elkies4q | k = 428843415, cover c1-5 | 10¹⁰ | 8 | 58 min 40 s | no point (the generator lies at max\|u\| ≈ 2·10¹¹) |
| elkies4 | same | 10¹⁰ | 8 | about 8 h | no point |
| elkies4q | k = 428843415, cover c1-5, half [0.5, 1) | 2·10¹² | 64 | about 12 h to the point (not a completed run; the other half ran at the same time on a second core) | point found |

These numbers compare our two programs only. Timings of other open-source searchers on the same inputs are welcome.
