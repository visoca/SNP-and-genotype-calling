## SNP and genotype calling with ANGSD

ANGSD is a suite of programmes for the analysis of NGS data. Among other things, it can be used to call SNP and genotypes using a varied range of approaches (e.g. it incorporates several SNP calling algorithms). The common theme of the programme is the use of a probabilistic framework that allows accounting for NGS-data uncertainties, especially with low-to-medium depth datasets.  

Here you can see a script to call SNPs in a similar fashion to bcftools or GATK, calling SNPs from three scaffolds in parallel, using 2 cores for each task:

```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l mem=2G
#$ -pe smp 2
#$ -t 1-3
#$ -j y
#$ -o angsd.log
#$ -N angsd

I=$(($SGE_TASK_ID-1))

# Regions (=scaffolds) that will be analysed
# There must be as many as tasks are specified above with '#$ -t'
REGIONS=(Hmel201001 Hmel201003 Hmel201008)

# Path to the genome
GENOME=/fastdata/$USER/varcal/genome/Hmel2.fa

# Path to directory with BAM files
ALIDIR=/fastdata/$USER/varcal/alignments

# Path to output directory
OUTDIR=/fastdata/$USER/varcal/angsd

# Create output directory
mkdir -p $OUTDIR >& /dev/null

# Log for this region
LOG=$OUTDIR/angsd-${REGIONS[$I]}.log

# Redirect all the following screen outputs to log file
exec > $LOG 2>&1


hostname
date
echo "=============================================================================="

# Load genomics software repository
source /usr/local/extras/Genomics/.bashrc

# Get list of bam files on a temporary file
BAMLIST=$(mktemp)
ls $ALIDIR/*.bam > $BAMLIST

# Call SNPs and genotypes with ANGSD
# -----------------------------------------------------------
angsd \
-nThreads 2 \
-bam $BAMLIST \
-ref $GENOME \
-r ${REGIONS[$I]}: \
-out $OUTDIR/angsd-${REGIONS[$I]} \
-doVcf 1 \
-GL 1 \
-doMaf 3 \
-doMajorMinor 4 \
-SNP_pval 1e-3 \
-doPost 1 \
-doCounts 1 \
-doGeno 11

# Index VCF (vcf.gz) file
# (required for subsequent conversion)
bcftools index $OUTDIR/angsd-${REGIONS[$I]}.vcf.gz

# Convert VCF to BCF
bcftools view $OUTDIR/angsd-${REGIONS[$I]}.vcf.gz -O b > $OUTDIR/angsd-${REGIONS[$I]}.bcf

# index BCF file
bcftools index $OUTDIR/angsd-${REGIONS[$I]}.bcf
# -----------------------------------------------------------

# remove temporary file
rm $BAMLIST

echo "=============================================================================="
date
```
