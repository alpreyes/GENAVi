#!/bin/env bash
#
#$ -cwd
#$ -v PATH=path/to/STAR-2.5.1b:path/to/preseq/preseq_v2.0
#$ -pe smp 24
#$ -l mem_free=36G
#$ -S /bin/bash

set -e

# This script assumes hg38 and gencode v26 genome build - that can be changed
# This script assumes that your reads are gzip'd
# This script assumes that your reads are named something like:
# <blahblahblah>.R1.fastq.gz and <blahblahblah>.R2.fastq.gz
# feed this script your file names as full paths, not just the names
# something like `readlink -f <filename>` will give you the full path of the file name

INFILE1="$1"
INFILE2="$2"
OUTDIR="$3"
CORES="24"

mkdir -p $OUTDIR

/path/to/STAR-2.5.1b/bin/Linux_x86_64_static/STAR \
        --runThreadN ${CORES} \
        --genomeDir 
/path/to/ref_genome_hg38_gencodev26 \
        --readFilesIn "${INFILE1}" "${INFILE2}" \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --outFileNamePrefix ${OUTDIR}/$(basename ${INFILE1%.R1.fastq.gz}) \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMmode Full \
        --outSAMattributes All \
        --outSAMunmapped Within

