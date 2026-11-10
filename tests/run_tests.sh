#!/usr/bin/env bash
# Regression tests for the compiled searches. Each case: a reduced four-cover (Q2, Q3 coefficients) with a known
# primitive point u; the search at bound H must print exactly that point (up to sign). Run from the repository root.
#   tests/run_tests.sh            # builds src/elkies4q (needs gcc with __float128/libquadmath) and src/elkies4
set -euo pipefail
cd "$(dirname "$0")/.."
gcc -O2 -Wall -Werror -o src/elkies4 src/elkies4.c -lm
HAVE_Q=1; gcc -O2 -Wall -Werror -o src/elkies4q src/elkies4q.c -lquadmath -lm 2>/dev/null || HAVE_Q=0
fail=0
check() {  # binary H VOL coef expected-point
  local bin=$1 H=$2 VOL=$3 coef=$4 want=$5
  local got; got=$("./src/$bin" "$H" "$VOL" 200 < "$coef" 2>/dev/null | awk '/^POINT/{print $2, $3, $4, $5}')
  local neg; neg=$(echo "$want" | awk '{print -$1, -$2, -$3, -$4}')
  if [ "$got" = "$want" ] || [ "$got" = "$neg" ]; then echo "ok   $bin H=$H $coef"; else echo "FAIL $bin H=$H $coef: got '$got', want '$want'"; fail=1; fi
}
while read -r coef H want; do
  [ -z "$coef" ] || [ "${coef:0:1}" = "#" ] && continue
  check elkies4 "$H" 8 "$coef" "$want"
  [ $HAVE_Q = 1 ] && check elkies4q "$H" 64 "$coef" "$want"
done < tests/cases.txt
[ $HAVE_Q = 1 ] || echo "note: elkies4q not built (no __float128 on this platform); only elkies4 tested"
exit $fail
