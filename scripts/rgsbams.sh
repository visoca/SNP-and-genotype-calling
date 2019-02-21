#!/bin/bash
#$ -l h_rt=2:00:00
#$ -l mem=2G
#$ -t 1-32
#$ -j y
#$ -o rgsbams.log
#$ -N rgsbams

I=$(($SGE_TASK_ID-1))

INDIR=/fastdata/$USER/alignments

ALIS=($(ls $INDIR/*.bam))

OUTDIR=/fastdata/$USER/alignments_rgs
mkdir -p $OUTDIR >& /dev/null

ALIIN=${ALIS[$I]}
ALIOUT=$OUTDIR/$(basename $ALIIN)

# get filename without path nor extension
ID=$(basename $ALIIN)
ID=${ID%.bam}

LOG=$OUTDIR/$ID.log

# Redirect all the following screen outputs to log file
exec > $LOG 2>&1

hostname
date
echo "=============================================================================="

# Load genomics software repository
source /usr/local/extras/Genomics/.bashrc 

cd $INDIR

picard AddOrReplaceReadGroups INPUT=$ALIIN SORT_ORDER=coordinate \
RGID=$ID RGLB=$ID RGPL=ILLUMINA RGSM=$ID RGPU=$ID CREATE_INDEX=true \
OUTPUT=$ALIOUT

samtools index $ALIOUT

echo "=============================================================================="
date
