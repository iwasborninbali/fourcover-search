#!/usr/bin/env bash
# End-to-end search for a rational point of E above chosen 2-covers.
#   1. search4.py: relative four-descent above each 2-cover, reduced quadrics for every soluble four-cover;
#   2. elkies4q (or elkies4): lattice search on Q2 = Q3 = 0 for each H in H_LIST, smallest first;
#   3. dense_solve.py: exact solution of any DENSE slice the search reports;
#   4. search4.py again: exact verification of every POINT, mapped down to E.
# Stops at the first H that yields a point. Completeness of the search is NOT certified (see README).
#
# Usage:  pipeline/run_world.sh WORLD_JSON "COVERS" "H_LIST" [VOL]
#   WORLD_JSON   written by pipeline/export_world.gp
#   COVERS       2-cover labels, e.g. "c1" or "c1 c2"
#   H_LIST       search bounds tried in order, e.g. "1e6 1e8 1e10"
#   VOL          box volume target for elkies4q (default 64; elkies4 uses its own default when called as ELKIES)
# Environment:
#   SAGE_PYTHON           Sage's Python (default: "sage -python")
#   ELLIPTIC_DESCENT_SRC  src/ directory of https://github.com/long-algorithms/elliptic-descent (commit 882a0dac)
#   ELKIES                search binary (default: src/elkies4q of this repository)
#   OUT                   output directory (default: runs/<world name>)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(dirname "$HERE")
W=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); COVERS=$2; HLIST=$3; VOL=${4:-64}
SAGE_PY=${SAGE_PYTHON:-"sage -python"}
ELKIES=${ELKIES:-$ROOT/src/elkies4q}
NAME=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["world"])' "$W")
OUT=${OUT:-runs/$NAME}; mkdir -p "$OUT"; export OUT_DIR=$OUT

for cv in $COVERS; do
  log=$OUT/descent-$cv.log
  [ -s "$log" ] || WORLD_JSON=$W $SAGE_PY "$HERE/search4.py" "$cv" 0 > "$log"
  python3 - "$log" "$OUT" "$cv" <<'EOF'
import json, sys
log, out, cv = sys.argv[1:]
for line in open(log):
    if line.startswith('{'):
        r = json.loads(line)
        if 'Q2' in r:
            open(f"{out}/cover-{cv}-{r['label']}.coef", 'w').write(' '.join(map(str, r['Q2'] + r['Q3'])) + '\n')
            print('soluble four-cover', cv, r['label'], 'xi_bits', r['xi_bits'], 'coeff_digits', r['reduced_coeff_digits'])
EOF
done

for H in $HLIST; do
  for cv in $COVERS; do
    for f in "$OUT"/cover-"$cv"-*.coef; do
      [ -e "$f" ] || continue
      lab=${f##*/cover-$cv-}; lab=${lab%.coef}; o=$OUT/search-$cv-$lab-H$H
      "$ELKIES" "$H" "$VOL" 200 < "$f" > "$o.out" 2> "$o.err"
      if grep -q ^DENSE "$o.out"; then $SAGE_PY "$HERE/dense_solve.py" "$f" "$o.out" > "$o.dense"; else : > "$o.dense"; fi
      { grep -h ^POINT "$o.out" "$o.dense" || true; } | awk -v l="$lab" '{print l, $2, $3, $4, $5, 1}' > "$o.points"
      if [ -s "$o.points" ]; then
        WORLD_JSON=$W $SAGE_PY "$HERE/search4.py" "$cv" 0 "$o.points" > "$o.verify"
        echo "FOUND cover $cv four-cover $lab H $H"
        python3 - "$o.verify" <<'EOF'
import json, sys
for line in open(sys.argv[1]):
    if line.startswith('{'):
        for v in json.loads(line).get('verified', []):
            if v.get('on_E'):
                print('POINT_ON_E', v['Q'][0], v['Q'][1], 'height', v['height_Q'])
EOF
        exit 0
      fi
      echo "NULL cover $cv four-cover $lab H $H ($(grep ^DONE "$o.out"))"
    done
  done
done
echo EXHAUSTED
