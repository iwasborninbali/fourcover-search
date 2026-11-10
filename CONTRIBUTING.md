# Contributing

Pull requests from people and from AI agents are welcome. If you are an agent, say so in the pull request and name who
operates you. We review agent contributions by the same rules as everyone else's.

## What helps most

- **Speed:** a faster `elkies4q` at the same bound H. Measure it on the benchmark in `bench/`.
- **Completeness:** an argument or implementation change that makes "every point with |u|₂ ≤ H" certifiable. Examples
  are rigorous curvature bounds, interval arithmetic for the boxes, and guaranteed discovery of all real components.
- **Portability:** a version for platforms without `__float128` (e.g. double-double arithmetic on arm64) that matches
  `elkies4q` on the tests.
- **Comparisons:** measured runs of other open-source point searchers on the same four-covers.
- **New results:** more curves solved with the pipeline. Add a folder under `results/` with the world file, the reduced
  four-cover, the point, and stdlib or PARI checks that someone else can run.

## Rules for a pull request

1. `tests/run_tests.sh` passes. CI runs it on Linux x86_64 with `REQUIRE_QUAD=1`, so `elkies4q` must build there. If you
   change search behaviour, add a test case to `tests/cases.txt`: a reduced four-cover, a bound H and a point that must
   be found. These are recovery tests, not completeness tests.
2. If the change is meant to speed up the search, include `bench/bench.sh` timings before and after, on the same machine.
3. Keep the search's fixed-width integer checks and successful `search4.py` map checks. Independently verify every released point with arbitrary-precision integer/rational arithmetic on the reduced quadrics and on E, and check that integer intermediates cannot overflow in the reported runs. A `FOUND` line alone is not successful verification.
4. Don't overstate. If a change rests on a heuristic, say so in the code comment and in the pull request.
5. Code is contributed under MIT; data and text under CC BY 4.0.

## Releases

Maintainers merge what passes review and tag a release (`vX.Y`). Each release is archived on Zenodo with a DOI, and
contributors are listed in the release notes and in `CITATION.cff`.
