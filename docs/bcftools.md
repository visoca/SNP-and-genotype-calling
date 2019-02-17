## 2. SNP and genotype calling with BCFtools
BCFtools is a very popular programme to call SNPs and genotypes (and also to manipulate and filter vcf/bcf files as we will see afterwards). SNP calling is a relatively intensive process, to speed things up we will be restricting variant calling to 3 scaffolds. Before calling SNPs, we have to decompress and index the genome using `samtools faidx`:
```bash
gzip -d genome/Hmel2.fa.gz
samtools faidx genome/Hmel2.fa
```
This will produce a .fai index file:
```bash
 ls -lh genome
```
>``-rw-r--r-- 1 myuser cs 76M Feb 13 17:11 Hmel2.fa``<br>
>``-rw-r--r-- 1 myuser cs 31K Feb 13 17:11 Hmel2.fa.fai``<br>

We are going to prepare now a batch script to call SNPs and genotypes for all 32 individuals. To speed things up, we will be using an [SGE array job](http://docs.hpc.shef.ac.uk/en/latest/parallel/JobArray.html) to call SNPs for three scaffolds (Hmel201001, Hmel201002, and Hmel201003) in parallel:

```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l mem=2G
#$ -t 1-3
#$ -m bea
#$ -M myemail@mail.com
#$ -j y
#$ -o bcftools.log
#$ -N bcftools

# internal SGE task index
I=$(($SGE_TASK_ID-1))

# Regions (=scaffolds) that will be analysed
# There must be as many as tasks are specified above with '#$ -t'
REGIONS=(Hmel201001 Hmel201003 Hmel201008)

# Path to the genome
GENOME=/fastdata/$USER/varcal/genome/Hmel2.fa

# Path to directory with BAM files
ALIDIR=/fastdata/$USER/varcal/alignments

# Path to output directory
OUTDIR=/fastdata/$USER/varcal/bcftools

# Create output directory
mkdir -p $OUTDIR >& /dev/null

# Log for this region
LOG=$OUTDIR/bcftools-${REGIONS[$I]}.log

# Redirect all the following screen outputs to log file
exec > $LOG 2>&1

hostname
date
echo "=============================================================================="

# Load genomics software repository
# (safer to load it here in case the worker nodes don't inherit the environment)
source /usr/local/extras/Genomics/.bashrc 

# Output BCF for this region
OUTBCF=$OUTDIR/bcftools-${REGIONS[$I]}.bcf
# Log for this region
LOG=$OUTDIR/bcftools-${REGIONS[$I]}.log

# Change to alignment dir
# This is necessary to avoid including the BAM filepaths the 
# samples ID in the BCF 
cd $ALIDIR

# Call SNPs and genotypes with bcftools
# -----------------------------------------------------------
#  mpileup
#   * -q 20: filter out alignments with mapping quality <20
#   * -Q 20: filter out bases with QS < 20
#   * -P ILLUMINA: use Illumina platform for indels
#   * -a FORMAT/DP,FORMAT/AD: output depth and allelic depth
#  call
#   * -m: use the multiallelic caller
#   * -v: output variants only
#   * -P 1e-6: prior on mutation rate: 
#   * -f GQ: ouput genotype quality
#   * -O b: output in BCF format

bcftools mpileup -Ou \
--max-depth 10000 \
-q 20 -Q 20 -P ILLUMINA \
-a FORMAT/DP,FORMAT/AD \
-f $GENOME \
-r ${REGIONS[$I]} *.bam | \
bcftools call -mv \
-P 1e-6 \
-f GQ \
-O b \
-o $OUTDIR/bcftools-${REGIONS[$I]}.bcf

# Index bcf file
bcftools index $OUTDIR/bcftools-${REGIONS[$I]}.bcf
# -----------------------------------------------------------

echo "=============================================================================="
date
```
Calling SNPs with bcftools is a two-step process. First, `bcftools mpileup` estimates genotype likelihoods at each genomic position with sequence data. Second, `bcftools call` identifies both variants and genotypes, i.e. makes the actual call. To avoid generating intermediate temporary files, the output of `bcftools mpileup` is *piped* to `bcftools call`. We are using a number of non-default options:
1. Estimating genotype likelihoods with `bcftools mpileup`:
   * `-q 20`: filter out sites with base quality (BQ) <20
   * `-Q 20`: filter out alignments with mapping quality (MQ) <20
   * `-P ILLUMINA`: specify the platform from which indels are collected (Illumina in this case)
   * `-a FORMAT/DP, FORMAT/AD`: output total depth (DP) and allelic depth (AD)
2. Calling SNPs, indels and genotypes with `bcftools call`:
   * `-m`: use the multiallelic caller, which implements a model that allows more than two alleles
   * `-v`: output only variants
   * `-P 1e-6`: use a relatively stringent mutation rate (the lower the more stringent)
   * `-f GQ`: output genotype quality (GQ)
   * `-O b`: output in BCF format, the binary version of VCF (more details below)

When you have finished editing the bash script, save it as `bcftools.sh`, make it executable with `chmod` and submit it to the job queue with `qsub`:
```bash
chmod +x bcftools.sh
qsub bcftools.sh
```
If all goes well, it should take no longer than a few minutes to get the jobs finished. You should now have the following files:
```ls -lh bcftools```

>``total 252K``<br>
>``-rw-r--r-- 1 myuser cs  95K Feb 17 01:57 bcftools-concat.bcf``<br>
>``-rw-r--r-- 1 myuser cs  287 Feb 17 01:57 bcftools-concat.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs  72K Feb 17 01:54 bcftools-Hmel201001.bcf``<br>
>``-rw-r--r-- 1 myuser cs  176 Feb 17 01:54 bcftools-Hmel201001.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs  505 Feb 17 01:54 bcftools-Hmel201001.log``<br>
>``-rw-r--r-- 1 myuser cs  26K Feb 17 01:54 bcftools-Hmel201003.bcf``<br>
>``-rw-r--r-- 1 myuser cs  158 Feb 17 01:54 bcftools-Hmel201003.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs  505 Feb 17 01:54 bcftools-Hmel201003.log``<br>
>``-rw-r--r-- 1 myuser cs  17K Feb 17 01:54 bcftools-Hmel201008.bcf``<br>
>``-rw-r--r-- 1 myuser cs  170 Feb 17 01:54 bcftools-Hmel201008.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs 9.0K Feb 17 01:55 bcftools-Hmel201008.log``<br>

and the content of the logfiles should look like this:

```less bcftools/bcftools-Hmel201001.log ```

>``sharc-node010.shef.ac.uk``<br>
>``Sun 17 Feb 01:54:29 GMT 2019``<br>
>``==============================================================================``<br>
>````<br>
>``  Your account is set up to use the Genomics Software Repository``<br>
>``     More info: http://soria-carrasco.staff.shef.ac.uk/softrepo``<br>
>````<br>
>``Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid``<br>
>``[mpileup] 32 samples in 32 input files``<br>
>``==============================================================================``<br>
>``Sun 17 Feb 01:54:41 GMT 2019``<br>

[Back to TOC](index.md)
