#!/bin/bash

# Raw, quick'n'dirty utility script to subset a big vcf file. Run this from grups/tests/test-data/vcf

function usage(){
    echo "Usage  : $0 <input-dir>"
    echo "Example: $0 ./g1k-phase3-v5b"
    exit
}
[[ $# -lt 1 ]] && usage;

IN_DIR="$1";
POP_REGEX="$2"
PANEL_DIR=${3:-$IN_DIR/integrated_call_samples_v3.20130502.ALL.panel}


POP_TAG="${POP_REGEX//|/-}"
OUT_DIR="${1%/}${POP_TAG:+-$POP_TAG}-filtered"

THREADS="$((`nproc`/2))"

echo INPUT     : $IN_DIR
echo OUTPUT DIR: $OUT_DIR
echo FILTER    : $POP_REGEX
echo PANEL     : $PANEL_DIR
echo THREADS   : $THREADS * 2
echo ""

mkdir -p ${OUT_DIR}
for i in ${IN_DIR}/*.vcf.gz; do 
    base="$(basename $i)";
    output="${OUT_DIR}/${base%.vcf.gz}${POP_TAG:+.$POP_TAG}.m2M2.snps.rmdup.vcf.gz";

    echo "Filtering: ${base} --> ${output}"

    # In case of duplicated, multiallelic SNP lines, this will keep the first instance, and exclude all others.
    #bcftools view --threads ${THREADS} -m2 -M2 --type snps $i | bcftools norm --threads ${THREADS} -Oz --rm-dup both -o ${output};

    # This on this other hand, will remove all the duplicated lines for these "false" biallelic SNP positions
    bcftools norm --threads ${THREADS} --multiallelics +snps $i \
	    | bcftools view --threads ${THREADS} \
	                    -m2 -M2 \
			    --type snps \
			    --samples-file <(cat "${PANEL_DIR}" | grep -P "${POP_REGEX}" | awk '{print $1}') \
			    -Oz -o ${output};

    tabix ${output};
done
