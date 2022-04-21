#!/bin/bash

BASEURL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/";

wget -qO- ${BASEURL} | grep -oP "ALL.chr[0-9]{1,2}.*vcf(.gz.tbi|.gz)" | sed "s/^/${BASEURL//\//\\/}/" > ./net/1000g-phase3-v5b-urls
