#!/bin/bash

set -eu -o pipefail

SCRIPT_DIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
TEST_DIR="${SCRIPT_DIR}/../test-data"
MIN_MQ=30
MIN_BQ=30


# ---- Horrid API
if [ $# -lt 2 ]; then 
    echo "$0 <URL_FILE> <TARGET_NAME> [THREADS]"
    echo ""
    echo "    - URL_FILE    (Required): simple, single column file containing the URLs leading to the bam required for fetching."
    echo "    - TARGET_NAME (Required): Which target file to use. Possible values:"
    echo "       --> HO    : ReichLab Humans Origin dataset"
    echo "       --> 1240K : ReichLab 1240K Target capture dataset"
    echo "    - THREADS     (Optional): Number of threads to use when fetching bams"
    echo "                              Default: `nproc`"
    echo ""
    echo "Example(s):"
    echo "    user@desktop:~$ $0 ./net/ART-test-data-urls.tsv 1240K"
    echo ""
    echo "      --> Expected output: `dirname $0`/test-data/pileup/ART-Bq30Q30-g1k_v37.pileup" 
    echo ""   
    exit
fi 

BAM_URLS="`pwd`/$1"
SNP_REGEX="$2"
    
# ---- Quick n' Dirty sanity checks:
[[ -f ${BAM_URLS} ]] || bail "${BAM_URLS} not found."
echo ${SNP_REGEX} | egrep -q '1240K|HO' || bail "Invalid TARGET_NAME: ${SNP_REGEX}" 

SNP_URL="https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V50/V50.0/SHARE/public.dir/v50.0_${SNP_REGEX}_public.snp"
REF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
FAI_URL="${REF_URL%.gz}.fai"

MAX_THREADS=`nproc`

function log(){
    # Basic logger, with 3 levels : INFO  - informational messages
    #                               WARN  - not fatal, but user should still get notified...
    #                               FATAL - Something went very wrong and should get fixed...
    local level message
    read -r level message <<< $(echo $1 $2) 
    local reset='\e[97m'
    declare -A levels=([INFO]="\e[32m" [WARN]="\e[33m" [FATAL]="\e[31m")
    local color=${levels["${level}"]}    
    echo -e "${color}[${level}]: ${message}${reset}"
}

function bail(){
    # Sends a log with FATAL level, accompanied with the line number, exit code and an optional message. 
    # Then, exit.. 
    local exit_code=$?
    local message=$1
    log FATAL "Line ${LINENO:-}: $message (exit code = ${exit_code})"
    exit $exit_code
}

function drawline(){
    # Draws a horizontal line of dashes, spanning the entire console length.
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function header(){
    # Calls drawline(), followed by a centered title
    local title=$1
    drawline
    title_len="${#title}"
    COLUMNS=$(tput cols)
    printf "%*s\n" "$(((title_len+COLUMNS)/2))" "$title"
}

# ----   I. Create test-data directory next to this script
header "Preparing workspace." 
BAM_DIR="${TEST_DIR}/bams"
REF_DIR="${TEST_DIR}/fasta"
SNP_DIR="${TEST_DIR}/targets"
PIL_DIR="${TEST_DIR}/pileups"

mkdir -p ${BAM_DIR} ${REF_DIR} ${SNP_DIR} ${PIL_DIR} || bail "Cannot create directories within ${TEST_DIR}" 

# ----  II. Download BAMS
header "Downloading Input BAM files"
pushd "${BAM_DIR}"
cat ${BAM_URLS} | grep -v "#" | xargs -n 1 -P ${MAX_THREADS} wget -Nq --show-progress || bail "Error when downloading bam"
popd

# ---- III. Download reference genome (fasta.gz + fasta.fai)
header "Downloading reference genome and index"
pushd "${REF_DIR}"
REF="${REF_DIR}/$(basename ${REF_URL%.gz})"
if [ ! -f "${REF}" ]; then 
    wget -q --show-progress -O- "${REF_URL}" | gunzip -cv > ${REF}
fi || bail "While downloading reference file: ${REF_URL}. Delete ${REF} and retry.\n"
wget -Nq --show-progress "${FAI_URL}"
popd 

# ----  IV. Download target sites
header "Downloading 1240K target sites"
pushd "${SNP_DIR}"
wget -Nq --show-progress "${SNP_URL}" 
popd

# ----   V. Perform mpileup
header "Running samtools mpileup"
pushd "${PIL_DIR}"
REGEX=$(basename ${BAM_URLS} | grep -oP '^.*?(?=[-_.])')
PILEUP="${PIL_DIR}/${REGEX}-Bq${MIN_MQ}Q${MIN_BQ}-${SNP_REGEX}-$(basename "${REF%.fasta}").pileup"
BAMS=$(cat ${BAM_URLS} | grep -v "#" | xargs basename -a | sed "s/^/${BAM_DIR//\//\\/}\//")
samtools mpileup -B -q${MIN_MQ} -Q${MIN_BQ}                                   \
	         -f ${REF}                                                    \
		 -l <(awk '{print $2, $4}' ${SNP_DIR}/$(basename ${SNP_URL})) \
		 ${BAMS}                                                      > ${PILEUP}&

# ----- Stupid progress bar.
FAI="${REF_DIR}/$(basename ${FAI_URL})"
sleep 2; echo ""
while [[ -n $(jobs -r) ]]; do 
    echo -e "\e[1A\e[KProgress: $(awk '
        FNR==NR{a[$1]=$2; next}
	($1 in a){printf "chr: %-3s - pos: %-10s - %.3f %", 
	                 $1,
			 a[$1],
			 (($1-1)/22+((a[$1]/$2)/22))*100}' <(tail -n-1 ${PILEUP}) ${FAI})";
    sleep 0.1;
done
wait
header "Done! "
log INFO "Output pileup is: ${PILEUP}"
exit
