#!/bin/bash

# Raw, quick'n'dirty utility script to subset a big vcf file. Run this from grups/tests/test-data/vcf

function usage(){
    echo "Usage  : $0 <nSNP>"
    echo "Example: $0 500"
    exit
}
[[ $# -lt 1 ]] && usage;

SNP=$1;

OUT_DIR="./g1k-phase3-v5b-first-${SNP}"
mkdir -p ${OUT_DIR}

for i in ./g1k-phase3-v5b/*.vcf.gz; do 
    base="$(basename $i)";
    output="${OUT_DIR}/${base%.vcf.gz}.first-${SNP}.vcf";
    bcftools view -h $i > $output ;
    bcftools view -H $i | head -n ${SNP}>> $output;
    bgzip $output;
    tabix ${output}.gz;
done
