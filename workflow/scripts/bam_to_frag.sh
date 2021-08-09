#!/bin/bash

FINAL_BAM_FILE=$1

FINAL_FRAG_FILE=$2
TMP_FRAG_FILE=$3

# --min_mapq 0 : no MAPQ constraints, since filtering is already done and everything should be accepted
# --max_distance 2000: This matches -X2000 argument for bowtie2, so should not filter any further
# --min_distance 10: minimum fragment length 
# these thresholds should not perform any further filtering on BAM-- need to verify

sinto fragments -b ${FINAL_BAM_FILE} -f ${TMP_FRAG_FILE} --min_mapq 0 --max_distance 2000 --min_distance 10 --barcodetag CB --nproc 25

cat ${TMP_FRAG_FILE} | sort -k1,1 -k2,2n | bgzip > ${FINAL_FRAG_FILE}

tabix -p bed ${FINAL_FRAG_FILE}

# rm ${TMP_FRAG_FILE}


