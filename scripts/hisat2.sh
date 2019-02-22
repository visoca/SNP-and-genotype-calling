#!/bin/bash
#$ -l h_rt=7:00:00
#$ -l mem=2G
#$ -pe smp 4
#$ -t 1-32
#$ -j y
#$ -o hisat2.log
#$ -N hisat2

I=$(($SGE_TASK_ID-1))

FQDIR=/usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/.original/butterflyRNAseq/trimmed_reads
OUTDIR=/fastdata/$USER/varcal/alignments_hisat2

mkdir -p $OUTDIR

READS1=($(find $FQDIR -name "*_1.fastq.gz"))
READS2=($(find $FQDIR -name "*_2.fastq.gz"))

ID=$(basename ${READS1[$I]})
ID=${ID%.trimA_1.fastq.gz}

ALI=$OUTDIR/$ID.bam
LOG=$OUTDIR/$ID.log

# Redirect all the following screen outputs to log file
exec > $LOG 2>&1

hostname
date
echo "=============================================================================="

# Load genomics software repository
source /usr/local/extras/Genomics/.bashrc 

GENOME=/fastdata/$USER/varcal/genome/Hmel2.fa

# RGS
RGPU=$(gzip -dc ${READS1[$I]} | head -n1 | awk -F: '{print $3"."$4}')
RGLB="Hmel"
RGID=$(echo $RGPU | perl -pe 's/ANX.*\./\./g')
RGPL="illumina"
RGSM=$ID

hisat2 \
--rg ID:$RGID --rg SM:$RGSM --rg PL:$RGPL --rg LB:$RGLB --rg PU:$RGPU \
-p 4 \ 
--rna-strandness RF \
-x ${GENOME%.fa} \
-1 ${READS1[$I]} -2 ${READS2[$I]} | \ 
samtools sort -m 6G | samtools view -b > $ALI 2> $LOG

samtools index $ALI

echo "=============================================================================="
date
