#!/bin/bash

FASTQ1_PATH=$1
FASTQ2_PATH=$2
BWT_IDX_PATH=$3

MULTIMAPPING=$4

BAM_RAW_PATH=$5
BAM_NO_MITO_PATH=$6
BAM_MITO_PATH=$7
FLAGSTAT_QC_PATH=$8

LOG_BWT_PATH=$9

THREADS=${10}

if [[ $MULTIMAPPING == 0 ]]
then
    bowtie2 -X2000 --mm -x $BWT_IDX_PATH --threads $THREADS -1 $FASTQ1_PATH -2 $FASTQ2_PATH --sam-append-comment \ 
    2> $LOG_BWT_PATH \
    | samtools view -u /dev/stdin \
    | samtools sort -o $BAM_RAW_PATH - 
else
    bowtie2 -k $((MULTIMAPPING+1)) -X2000 --mm -x $BWT_IDX_PATH --threads $THREADS -1 $FASTQ1_PATH -2 $FASTQ2_PATH --sam-append-comment \ 
    # 2> $LOG_BWT_PATH \
    # | samtools view -u /dev/stdin \
    # | samtools sort -o $BAM_RAW_PATH - 
fi

printf "asdfasdf\n"

# samtools sort -n --@ $THREADS -O SAM $BAM_RAW_PATH | SAMstats --sorted_sam_file - --outf $FLAGSTAT_QC_PATH

# printf "ddddd\n"

# samtools idxstats $BAM_RAW_PATH | cut -f 1 | grep -v -P "^chrM$" | xargs samtools view $BAM_RAW_PATH -@ $THREADS -b> $BAM_NO_MITO_PATH
# samtools idxstats $BAM_RAW_PATH | cut -f 1 | grep -P "^chrM$" | xargs samtools view $BAM_RAW_PATH -@ $THREADS -b> $BAM_MITO_PATH

