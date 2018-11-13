#!/bin/env bash
#
#$ -cwd
#$ -v PATH=path/to/preseq/preseq_v2.0
#$ -pe smp 4
#$ -l mem_free=12G
#$ -S /bin/bash

set -e


if [[ $1 =~ "clean" ]] ; then
        EXTEN="clean${1#*.clean}"
else
        EXTEN=${1:0-3}
fi

if [[ $EXTEN =~ "sam" ]] || [[ $EXTEN =~ "bam" ]] ; then
        mkdir -p $(dirname $1)/metrics
        mkdir -p $(dirname $1)/intermediate
        mkdir -p $(dirname $1)/original

        STORE="$(dirname $1)/intermediate"
        FILE=$(basename $1)
        FILESTUB=${FILE%.${EXTEN}}
        echo $FILESTUB starting

        if [[ ! -s ${FILESTUB}.clean.bam ]] && [[ ! -s ${STORE}/${FILESTUB}.clean.bam ]] ; then
                echo "executing picardtools CleanSam"
                /path/to/picardtools CleanSam \
                        TMP_DIR=/path/to/scratch/directory \
                        VERBOSITY=WARNING \
                        INPUT=$1 \
                        OUTPUT=${FILESTUB}.clean.bam &&
                mv $1 $(dirname $1)/original/
        fi
        if [[ ! -s ${FILESTUB}.clean.sorted.bam ]] && [[ ! -s ${STORE}/${FILESTUB}.clean.sorted.bam ]] ; then
                echo "executing picardtools SortSam"
                /path/to/picardtools SortSam \
                        TMP_DIR=/path/to/scratch/directory \
                        VERBOSITY=WARNING \
                        COMPRESSION_LEVEL=0 \
                        INPUT=${FILESTUB}.clean.bam \
                        OUTPUT=${FILESTUB}.clean.sorted.bam \
                        CREATE_INDEX=true \
                        SORT_ORDER=coordinate &&
                mv ${FILESTUB}.clean.bam ${STORE}/ &&
                mv ${FILESTUB}.clean.sorted.bai ${STORE}/${FILESTUB}.clean.sorted.bam.bai
        fi
        if [[ ! -s ${FILESTUB}.clean.sorted.mdups.bam ]] && [[ ! -s ${STORE}/${FILESTUB}.clean.sorted.mdups.bam ]] ; then
                echo "executing picardtools MarkDuplicates"
                /path/to/picardtools MarkDuplicates \
                        TMP_DIR=/path/to/scratch/directory \
                        VERBOSITY=WARNING \
                        INPUT=${FILESTUB}.clean.sorted.bam \
                        OUTPUT=${FILESTUB}.clean.sorted.mdups.bam \
                        CREATE_INDEX=true \
                        ASSUME_SORTED=true \
                        METRICS_FILE=$(dirname $1)/metrics/${FILESTUB} &&
                mv ${FILESTUB}.clean.sorted.bam ${STORE}/ &&
                mv ${FILESTUB}.clean.sorted.mdups.bai ${FILESTUB}.clean.sorted.mdups.bam.bai
        fi

        FINALBAM=${FILESTUB}.clean.sorted.mdups.bam
        DUPSBAM="${STORE}/${FILESTUB}.clean.sorted.bam"
        PREFIX="$(dirname $1)/metrics/${FILESTUB}"
        BED="/path/to/gencode.v26.annotation.bed"
        source activate rna-seqc
        echo 'executing: samtools mapped reads only' ; samtools view -bh -f 3 -F 4 $FINALBAM > ${FINALBAM%.bam}.mapped.bam &
        export samtools_PID=$!
        echo "samtools view has pid: $samtools_PID"
        echo 'calculating metrics now'
        echo 'executing: geneBody_coverage.py' ; geneBody_coverage.py -i $FINALBAM -r $BED -o $PREFIX &
        export geneBody_PID=$!
        echo "geneBody_coverage.py has pid: $geneBody_PID"
        echo 'executing: tin.py' ; tin.py -i $FINALBAM -r $BED -s > ${PREFIX}.tin &
        export tin_PID=$!
        echo "tin.py has pid: $tin_PID"
        echo 'executing: inner_distance.py' ; inner_distance.py -i $FINALBAM -r $BED -o $PREFIX
        echo 'executing: junction_annotation.py' ; junction_annotation.py -i $FINALBAM -r $BED -o $PREFIX > ${PREFIX}.junction_annotation.txt
        echo 'executing: junction_saturation.py' ; junction_saturation.py -i $FINALBAM -r $BED -o $PREFIX
        echo 'executing: read_distribution.py' ; read_distribution.py -i $DUPSBAM -r $BED > ${PREFIX}.read_distribution.txt
        echo 'executing: read_duplication.py' ; read_duplication.py -i $FINALBAM -o $PREFIX
        echo 'executing: read_GC.py' ; read_GC.py -i $FINALBAM -o $PREFIX
        echo 'executing: picard CollectRnaSeqMetrics' ; picardtools CollectRnaSeqMetrics $(printf "I=%s O=%s " $FINALBAM "${PREFIX}.picardRNA.metrics") \
                REF_FLAT=/path/to/gencode.v26.annotation.refFlat \
                VERBOSITY=WARNING \
                STRAND=NONE TMP_DIR=/path/to/scratch/directory \
                RIBOSOMAL_INTERVALS=/path/to/hg38.rRNA.interval_list
        wait $samtools_PID
        echo 'executing: preseq c_curve' ;
        preseq c_curve -l 5000000 -B -P -o ${PREFIX}.c_curve ${FINALBAM%.bam}.mapped.bam
        echo 'executing: preseq lc_extrap' ;
        preseq lc_extrap -e 100000000 -l 5000000 -B -P -o ${PREFIX}.lc_extrap ${FINALBAM%.bam}.mapped.bam
        wait $geneBody_PID
        wait $tin_PID
        echo 'all done'
        source deactivate
else
        echo "file does not end in either .sam or .bam"
fi

