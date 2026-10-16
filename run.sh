#!/bin/sh
# usage: ./run.sh orbit_count|line_count FIRST LAST [threads]
# Appends "n a(n) seconds" lines to results/<program>.txt.
prog=$1; first=$2; last=$3; threads=${4:-3}
for n in $(seq "$first" "$last"); do
    start=$(date +%s)
    out=$(./"$prog" "$n" "$threads") || { echo "$prog $n failed" >&2; exit 1; }
    echo "$out $(( $(date +%s) - start ))" >> "results/$prog.txt"
done
