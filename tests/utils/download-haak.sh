#!/bin/bash

SCRIPT_DIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))

HAAK_TSV="${SCRIPT_DIR}/net/Haak-test-data-urls.tsv"

BQ=30
MQ=30
REFERENCE="../../fasta/human_g1k_v37.fasta"
HAAK_TARGETS="../../targets/Haak-2015_390K.snp"

for line in $(awk 'BEGIN{OFS=","}{print $9, $10}' ${HAAK_TSV} ); do 
    IFS=',';
    set -- $line
    echo $1 $2
    output_bam=$(basename $(echo $1 | sed -E "s/S0([0-9]{3})/\\$2/"))
    if [ ! -f ${output_bam} ]; then 
        wget $1 -O $output_bam
    fi
done

#conda create -n mapdamage-2.0.6 -c bioconda -c r mapdamage2=2.0.6 r-base>=3.5
#install.packages(c('inline', 'ggplot2', 'gam', 'Rcpp', 'RcppGSL'))
# mamba install -c conda-forge r-rcppgsl
#  mamba install -c conda-forge r-ggplot2
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mapdamage-2.0.6

MAX_PROCS=`nproc`
PROCS=0

for bam in ./*.bam ; do
    ((PROCS++))
    echo "$bam | n-thread: $PROCS"
    mapDamage -i $bam -r ../../fasta/human_g1k_v37.fasta --rescale &
    if ((PROCS%MAX_PROCS==0)); then
        wait
        PROCS=0
    fi
done

for bam in results_I01*/*.bam; do
    echo $bam;
done > Esperstedt.bamset

samtools mpileup -q $BQ -Q $MQ -f ${REFERENCE} -l <(awk '{print $2, $4}' ${HAAK_TARGETS}) -b Esperstedt.bamset > Esperstedt.q30Q30.390K.pileup

# --------------------------------------
for bam in results_I04*/*.bam; do
    echo $bam;
done > Els-Troc.bamset

samtools mpileup -q $BQ -Q $MQ -f ${REFERENCE} -l <(awk '{print $2, $4}' ${HAAK_TARGETS}) -b Els-Troc.bamset > Els-Troc.q30Q30.390K.pileup

