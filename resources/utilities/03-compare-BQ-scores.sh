#!/usr/bin/env bash

BAM_COL=11
PYTHON_CMD=$(cat << EOM
import sys
from math import floor
PHRED_SHIFT=33
for line in sys.stdin:
    line    = line.strip('\n').split(" ")
    left    = line[0]
    right   = line[1]
    cmp_len = min(len(left), len(right))
    print("-------- Compoarison len: ", cmp_len)
    for (i, (a, b)) in enumerate(zip(left, right)):
        if a != b:
            relpos = floor(cmp_len/2) - floor(abs(i - (cmp_len/ 2)))
            print(a, ord(a)-PHRED_SHIFT, b, ord(b)-PHRED_SHIFT, relpos)
EOM
)

paste <(samtools view --no-header $1 |awk -v col=$BAM_COL '{print $col}') <(samtools view --no-header $2 | awk -vcol=$BAM_COL '{print $col}') | awk '{OFS=" "}($1!=$2){print $1, $2}' | python -c "${PYTHON_CMD}" | awk '{before_scores+=$2; after_scores+=$4; avg_relpos+=$5 }END{print "Average BQ-score before rescaling:", before_scores/NR; print "Average BQ-score after rescaling :", after_scores/NR; print "Average relative position:", avg_relpos/NR}'
