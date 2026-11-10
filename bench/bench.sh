#!/usr/bin/env bash
# Benchmark: re-find the generator's four-cover point for y^2 = x^3 + 730033053750 at H = 3e9 (about 20 min on one core).
# Usage: bench/bench.sh [binary] [VOL]     e.g. bench/bench.sh src/elkies4q 64
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=${1:-src/elkies4q}; VOL=${2:-64}
WANT="2045913748 -1557028403 -858328586 334741533"
start=$(date +%s)
got=$("$BIN" 3e9 "$VOL" 200 < results/k730033053750/cover-c1-156.coef 2>/dev/null | tee /dev/stderr | awk '/^POINT/{print $2, $3, $4, $5}')
end=$(date +%s)
neg=$(echo "$WANT" | awk '{print -$1, -$2, -$3, -$4}')
if [ "$got" = "$WANT" ] || [ "$got" = "$neg" ]; then ok=found; else ok="NOT FOUND (got '$got')"; fi
echo "bench k=730033053750 H=3e9 VOL=$VOL $BIN: $ok, wall $((end - start)) s on $(uname -m)"
