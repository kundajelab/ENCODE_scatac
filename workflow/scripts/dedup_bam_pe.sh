#!/bin/bash

RAW_BAM_FILE=$1

FINAL_BAM_FILE=$2
FINAL_BAM_INDEX_FILE=$3
FINAL_BAM_FILE_MAPSTATS=$4
DUP_FILE_QC=$5
PBC_FILE_QC=$6

TEMP_PREFIX=$7

MULTIMAPPING=$8
MMP_PATH=$9

THREADS=${10}

# =============================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# Only keep properly paired reads
# Obtain name sorted BAM file
# ==================
FILT_BAM_PREFIX="${TEMP_PREFIX}.filt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_PREFIX="${FILT_BAM_PREFIX}.nmsrt"
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
TMP_FILT_FIXMATE_BAM_FILE="${TMP_FILT_BAM_PREFIX}.fixmate.bam"

if [[ $MULTIMAPPING == 0 ]]
then
    MAPQ_THRESH=30
    samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | samtools sort -n -o ${TMP_FILT_BAM_FILE} -
    samtools fixmate -r ${TMP_FILT_BAM_FILE} ${TMP_FILT_FIXMATE_BAM_FILE}
else
    samtools view -F 524 -f 2 -u ${RAW_BAM_FILE} | samtools sort -n -o ${TMP_FILT_BAM_FILE} -
    samtools view -h ${TMP_FILT_BAM_FILE} | python ${MMP_PATH} -k $MULTIMAPPING --paired-end | samtools fixmate -r /dev/stdin ${TMP_FILT_FIXMATE_BAM_FILE}
fi

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM


samtools view -F 1804 -f 2 -u ${TMP_FILT_FIXMATE_BAM_FILE} | samtools sort -o ${FILT_BAM_FILE} -

rm ${TMP_FILT_FIXMATE_BAM_FILE}

rm ${TMP_FILT_BAM_FILE}

# =============
# Mark duplicates
# =============
# module add picard-tools/1.126

TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
MARKDUP="/srv/gs1/software/picard-tools/1.126/MarkDuplicates.jar"
# DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"

java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false BARCODE_TAG=CB

mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

# ============================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ============================
# FINAL_BAM_PREFIX="${OFPREFIX}.nodup"
# FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
# FINAL_BAM_INDEX_FILE="${FINAL_BAM_FILE}.bai"
# FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}


# Index Final BAM file
samtools index ${FINAL_BAM_FILE}

samtools sort -n -@ ${THREADS} ${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${FINAL_BAM_FILE_MAPSTATS}

# =============================
# Compute library complexity
# =============================
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statistics

# module add bedtools/2.26

# PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair


samtools sort -n ${FILT_BAM_FILE} -o ${TEMP_PREFIX}.srt.tmp.bam
bedtools bamtobed -bedpe -i ${TEMP_PREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
rm ${OFPREFIX}.srt.tmp.bam


rm ${FILT_BAM_FILE}
