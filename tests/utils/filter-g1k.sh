#!/bin/bash

# Raw, quick'n'dirty utility script to subset a big vcf file. Run this from grups/tests/test-data/vcf

function usage(){
    echo "Usage  : $0 <input-dir>"
    echo "Example: $0 ./g1k-phase3-v5b"
    exit
}
[[ $# -lt 1 ]] && usage;

IN_DIR="$1";
OUT_DIR="${1%/}-filtered"

THREADS="$((`nproc`/2))"
echo $THREADS
mkdir -p ${OUT_DIR}

for i in ${IN_DIR}/*.vcf.gz; do 
    base="$(basename $i)";
    output="${OUT_DIR}/${base%.vcf.gz}.m2M2.snps.rmdup.vcf.gz";

    echo "Filtering: ${base} --> ${output}"
    bcftools view --threads ${THREADS} -m2 -M2 --type snps $i | bcftools norm --threads ${THREADS} -Oz --rm-dup both -o ${output};
    tabix ${output};
done
