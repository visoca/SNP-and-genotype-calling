## SNP and genotype calling with ANGSD

ANGSD is a suite of programmes for the analysis of NGS data. Among other things, it can be used to call SNP and genotypes using a varied range of approaches (e.g. it incorporates several SNP calling algorithms). The common theme of the programme is the use of a probabilistic framework that allows accounting for NGS-data uncertainties, especially with low-to-medium depth datasets. It is less polished than bcftools and GATK and the documentation is not as good. However, it enables a number of other analyses, such as the estimation of site frequency spectra, ABBABABA tests, the estimation of several population genetics statistics such as FST, Thetas, Tajima's D, etc.



Here you can see a script to call SNPs in a similar fashion to bcftools or GATK, calling SNPs from three scaffolds in parallel, using 2 cores for each task. As with the previous exercises, you will need to have the indexed reference genome and a file with the BAM files.

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
-minMapQ 20 \
-minQ 20 \
-doVcf 1 \
-GL 1 \
-doMaf 3 \
-doMajorMinor 4 \
-SNP_pval 1e-6 \
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
The options used are:
* `-nThreads 2`: Use two cores for computations
* `-bam $BAMLIST`: List of BAM files (one per sample)
* `-ref $GENOME`: Reference sequence
* `-r ${REGIONS[$I]}:`: Regions to call variants from (notice
* `-out $OUTDIR/angsd-${REGIONS[$I]}` output base name
* `-minQ 20`: filter out bases with quality (BQ) <20
* `-minMapQ 20`: filter out alignments with mapping quality (MQ) <20
* `-doVcf 1`: output in VCF format
* `-GL 1`: estimate genotype likelihoods using bcftools method (2-GATK, 3-SOAPsnp)
* `-doMaf 3`: estimate allele frequencies using both
* `-doMajorMinor 4`: Use reference allele as major (similar to bcftools and GATK)
* `-SNP_pval 1e-6`: p-value for calling SNPs
* `-doPost 1`: estimate genotype posterior probabilities using allele frequencies as priors (2-use uniform prior)
* `-doCounts 1`: do and output allele counts
* `-doGeno 8`: call genotypes and write posterior probability of genotypes

You will see it generates several other files apart from the VCF/BCF files:
```bash
ls -lh angsd
```
>``total 688K``<br>
>``-rw-r--r-- 1 bo1vsx bo 8.9K Feb 18 05:07 angsd-Hmel201001.arg``<br>
>``-rw-r--r-- 1 bo1vsx bo 175K Feb 18 05:07 angsd-Hmel201001.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo  154 Feb 18 05:07 angsd-Hmel201001.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo  76K Feb 18 05:07 angsd-Hmel201001.geno.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo 1.7K Feb 18 05:07 angsd-Hmel201001.log``<br>
>``-rw-r--r-- 1 bo1vsx bo  16K Feb 18 05:07 angsd-Hmel201001.mafs.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo 177K Feb 18 05:07 angsd-Hmel201001.vcf.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo  179 Feb 18 05:07 angsd-Hmel201001.vcf.gz.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo 8.9K Feb 18 05:07 angsd-Hmel201003.arg``<br>
>``-rw-r--r-- 1 bo1vsx bo  39K Feb 18 05:07 angsd-Hmel201003.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo  138 Feb 18 05:07 angsd-Hmel201003.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo  22K Feb 18 05:07 angsd-Hmel201003.geno.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo 1.7K Feb 18 05:07 angsd-Hmel201003.log``<br>
>``-rw-r--r-- 1 bo1vsx bo 5.1K Feb 18 05:07 angsd-Hmel201003.mafs.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo  39K Feb 18 05:07 angsd-Hmel201003.vcf.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo  163 Feb 18 05:07 angsd-Hmel201003.vcf.gz.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo 8.9K Feb 18 05:07 angsd-Hmel201008.arg``<br>
>``-rw-r--r-- 1 bo1vsx bo  19K Feb 18 05:07 angsd-Hmel201008.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo  147 Feb 18 05:07 angsd-Hmel201008.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo 8.8K Feb 18 05:07 angsd-Hmel201008.geno.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo 1.6K Feb 18 05:07 angsd-Hmel201008.log``<br>
>``-rw-r--r-- 1 bo1vsx bo 2.5K Feb 18 05:07 angsd-Hmel201008.mafs.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo  19K Feb 18 05:07 angsd-Hmel201008.vcf.gz``<br>
>``-rw-r--r-- 1 bo1vsx bo  174 Feb 18 05:07 angsd-Hmel201008.vcf.gz.csi``<br>

In this case, `geno.gz` files contain genotype posterior probabilities (3 values per genotype) and `mafs.gz` contain allele frequencies. However, other files may be generated depending on the options used.

As you can see the content of the VCF/BCF file differs quite a bit from bcftools and GATK: only the scaffolds/chromosomes analysed are included in the header, there are no quality scores nor allele counts (AC and AN) for the SNPs, and only biallelic SNPs are called. Some of these fields can be added a posteriori with bcftools, for example for AC and AC:
```bash
bcftools +fill-tags angsd/angsd-Hmel201001.bcf -- -t AC,AN | less -S
```
Also you can see the number of called SNPs tend to be higher than with bcftools or GATK:
```bash
bcftools view -H angsd/angsd-Hmel201001.bcf | wc -l
```
Because of that, applying filters based on depths or fraction of samples with data is recommended. 
