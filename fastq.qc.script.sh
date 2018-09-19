#!/bin/env bash
#
#$ -v PATH=/path/to/FastQC
#$ -l mem_free=8G
#$ -S /bin/bash
#$ -cwd

set -e
PROGNAME=$(basename $0)

function error_exit
{

#   ----------------------------------------------------------------
#   Function for exit due to fatal program error
#       Accepts 1 argument:
#           string containing descriptive error message
#   ----------------------------------------------------------------


    echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
    exit 1
}

#source activate rna-seqc

FASTQ=$1
OUT_DIR=$2

/path/to/FastQC/fastqc -o $OUT_DIR -d /path/to/scratch/directory -f fastq $FASTQ
/path/to/fastq_screen_v0.6.3/fastq_screen --outdir $OUT_DIR $FASTQ

#source deactivate

