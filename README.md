**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# SNP and genotype calling
#### Victor Soria-Carrasco

The aim of this practical is to learn how to call single nucleotide polymorphism (SNPs) and genotypes, that is the process of identifying variable sites and determining the genotype for each individual at each site. We will be using a dataset of whole genome sequence data of 32 individuals of *Heliconius melpomene*. After calling SNPs, we will do some subsetting and filtering and will carry out a few example analyses.

## Table of contents
1. [SNP and genotype calling with BCFtools](#SNP-and-genotype-calling-with-BCFtools)
2. [VCF and BCF format](#VCF-and-BCF-format)
3. [SNP and genotype calling with GATK](#SNP-and-genotype-calling-with-GATK)
4. [Operations with BCF files](#Operations-with-BCF-files)

Extras:<br>
  * [SNP and genotype calling with ANGSD](#ANGSD-SNP-and-genotype-calling)<br>
  * [Population structure with NGSADMIX](#Population-structure-with-NGSADMIX)<br>
  * [PCA of genoypes with R](#PCA-of-genoypes-with-R)<br>
  * [Differentiation genome scans using pairwise FSTs](#Differentiation-genome-scans-using-pairwise-FSTs)<br>

---
---

## Initial set up
First of all, please remind that this tutorial must be run using an interactive session in ShARC. For that, you should log in into ShARC with `ssh`, and then request an interactive session with `qrsh`. Your shell prompt should show `sharc-nodeXXX` (XXX being a number between 001 and 172) and not `@sharc-login1` nor `@sharc-login2`.

For this particular tutorial, we are going to create and work on a directory called `varcal` in your /fastdata/$USER directory:
```bash
cd /fastdata/$USER
mkdir varcal
cd varcal
```
Remember you can check your current working directory anytime with the command `pwd`.
It should show something like:
```bash
pwd
```
>`´/fastdata/myuser/varcal`´<br>

## Programmes
There are many different tools for SNP calling, we are going to use three of the most popular: [bcftools](http://www.htslib.org/), [GATK](https://software.broadinstitute.org/gatk/), and, optionally,[ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). Additionally, we will use bcftools and awk for filtering files and picard tools for adding information about read groups to BAM files. All of them are already installed as part of the [Genomics Software Repository](http://soria-carrasco.staff.shef.ac.uk/softrepo/).

## Data
We are going to use the directory `/fastdata/$USER/varcal` as the base working directory for this tutorial, let's create and change to it:
```bash
mkdir -p /fastdata/$USER/varcal
cd /fastdata/$USER/varcal
```
Now we need to copy the alignments in BAM format produced in the previous sessions to `/fastdata/$USER/varcal/alignments`. If you don't have such data, there are files are available in ` /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/data/varcal/bams`, e.g.:
```bash
cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/data/varcal/alignments ./
```
It is very important that all the BAM files are indexed, which can be done the following way:
```bash
ls alignments/*.bam | xargs -I {} sh -c 'samtools index {}'
```
The files should look like this:
```bash
ls -lh alignments
```
>``-rw-r--r--   1 myuser bo 1.9G Feb 13 04:29 60A.bam``<br>
>``-rw-r--r--   1 myuser bo 482K Feb 13 04:30 60A.bam.bai``<br>
>``-rw-r--r--   1 myuser bo 1.7G Feb 13 04:29 60I.bam``<br>
>``-rw-r--r--   1 myuser bo 468K Feb 13 04:29 60I.bam.bai``<br>
>``-rw-r--r--   1 myuser bo 2.1G Feb 13 04:31 61A.bam``<br>
>``-rw-r--r--   1 myuser bo 586K Feb 13 04:31 61A.bam.bai``<br>

Lastly, we also need to download the reference genome, you can do it using [wget](https://www.gnu.org/software/wget/manual/wget.html) (remember that you can also read the manual on the command line with `man wget`):
```bash
mkdir genome
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2_-_scaffolds.fa.gz -O genome/Hmel2.fa.gz
```

---

## 1. SNP and genotype calling with BCFtools
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
REGIONS=(Hmel201001 Hmel201002 Hmel201003)

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

>``total 1.6M``<br>
>``-rw-r--r-- 1 myuser cs  72K Feb 15 18:46 bcftools-Hmel201001.bcf``<br>
>``-rw-r--r-- 1 myuser cs  176 Feb 15 18:46 bcftools-Hmel201001.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs  505 Feb 15 18:46 bcftools-Hmel201001.log``<br>
>``-rw-r--r-- 1 myuser cs 1.4M Feb 15 18:50 bcftools-Hmel201002.bcf``<br>
>``-rw-r--r-- 1 myuser cs 1.5K Feb 15 18:50 bcftools-Hmel201002.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs  505 Feb 15 18:50 bcftools-Hmel201002.log``<br>
>``-rw-r--r-- 1 myuser cs  26K Feb 15 18:46 bcftools-Hmel201003.bcf``<br>
>``-rw-r--r-- 1 myuser cs  158 Feb 15 18:46 bcftools-Hmel201003.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs  505 Feb 15 18:46 bcftools-Hmel201003.log``<br>

and the content of the logfiles should look like this:

```less bcftools/bcftools-Hmel201001.log ```

>``sharc-node114.shef.ac.uk``<br>
>``Fri 15 Feb 18:46:38 GMT 2019``<br>
>``==============================================================================``<br>
><br>
>``  Your account is set up to use the Genomics Software Repository``<br>
>``     More info: http://soria-carrasco.staff.shef.ac.uk/softrepo``<br>
><br>
>``Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid``<br>
>``[mpileup] 32 samples in 32 input files``<br>
>``==============================================================================``<br>
>``Fri 15 Feb 18:46:50 GMT 2019``<br>

## 2. VCF and BCF format
One of the most popular text-file formats for storing genetic variation information is the [Variant Call Format (VCF)](http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)). There are different versions of the format, but the core elements are the same. You can find a full description of the latest iteration [here](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf).Due to the large amount of data usually involved, files tend to be stored compressed. BCF is the binary version of VCF, which keeps the same information, but compressed and indexed. It  is designed to be much more efficient to process and store large amounts of data. The practical downside is that the contents can only be accessed with `bcftools view`.

The VCF format is composed of meta-information lines (prefixed wih ##), a header line (prefixed with #), and data lines each containing information about a position in the genome and genotype information on samples for each position (text fields separated by tabs). There are 8 fixed field lines, usually followed by an additional FORMAT field and an arbitrary number of sample fields when genotypes are included.

**Table 1** Description of VCF fields

Column Number| Title | Description
:--|:--|:--
1 | CHROM | Chromosome/scaffold of the reference genome
2 | POS | Position on the chromosome/scaffold given in CHROM (1-based)
3 | ID | ID for the SNP/INDEL (Here always '.')
4 | REF | Reference allele (usually the allele on reference genome)
5 | ALT | Alternate alleles (separated by commas);
6 | QUAL | Variant quality score ([phred-scale](https://gatkforums.broadinstitute.org/gatk/discussion/4260/how-should-i-interpret-phred-scaled-quality-scores))
7 | FILTER | Filter status ('.' if unfiltered or hard-filtered, PASS or a semicolon-separated list of codes for filters that fail if soft-filtered)
8 | INFO | List of variant annotations delimited by ';' and as specified in the meta-information lines
9 | FORMAT | Describes the format for the sample genotype information in column 10 onward
10 to end | Sample id (from SM tag in BAM file RG string) | Sample genotype information. Usually stores the genotype call (GT), genotype likelihoods (GL or PL (phred-scaled)), depth (DP), allelic depth (AD), and genotype quality (GQ).
 
Let's now have a look at our files:
```bash
bcftools view bcftools/Hmel201001.bcf | less -S
```
The file starts with a few lines specifying the VCF and bcftools versions, the command used with `bcftools mpileup`, followed by a long list of the >700 scaffolds that comprise the genome assembly (starting with `##contig`), the descriptions of the fields for ALT, INFO, and FORMAT, and any additional commands that were executed to produce this file (`bcftools call` in this case).

To focus only on the header of the SNPs, the flags `-h` and `-H` can be used:
```bash
bcftools view -h bcftools/Hmel201001.bcf | less -S
bcftools view -H bcftools/Hmel201001.bcf | less -S
```

## 3. SNP and genotype calling with GATK
GATK is another popular alternative. The algorithms used are more complex than those of bcftools, which makes the process of SNP calling slower. You can find how ``HaplotypeCaller`` - the caller we will be using in this practical - works [here](https://software.broadinstitute.org/gatk/documentation/article?id=4148). Another advange is its good documentation, with frequently updated guides on *Best Practices*. For example, you can find a document on the [GATK Best Practices for calling variants in RNASeq](https://software.broadinstitute.org/gatk/documentation/article.php?id=3891). However, we are not following all these recommendations here, because of time limitations.  

Before calling SNPs, we have to create the dictionary file for the genome:
```bash
gatk CreateSequenceDictionary -R genome/Hmel2.fa -O genome/Hmel2.dict
```
This will produce a dictionary file:
```bash
 ls -lh genome
```
>``-rw-r--r-- 1 myuser cs 92K Feb 13 19:50 Hmel2.dict``<br>
>``-rw-r--r-- 1 myuser cs 76M Feb 13 17:11 Hmel2.fa``<br>
>``-rw-r--r-- 1 myuser cs 31K Feb 13 17:11 Hmel2.fa.fai``<br>

__Important note:__ GATK is more picky than bcftools with regards to the content of the BAM files and requires to have read groups information appropriately encoded. The BAM files that will be use below have been produced to include such information, but it is something to bear in mind if plannig to use GATK. I have included in this repository two examples scripts on (1) how to run an aligner such as [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) to include read groups in the output BAM files ([hisat2.sh](https://github.com/visoca/SNP-and-genotype-calling/blob/master/scripts/hisat2.sh)) and (2) how to run [picard tools](https://broadinstitute.github.io/picard/) to add read groups to the BAM files *a posteriori* ([rgsbams.sh](https://github.com/visoca/SNP-and-genotype-calling/blob/master/scripts/rgsbams.sh)). 

Now, let's prepare a batch script to run GATK on ShARC. As before, we will use an SGE array to call SNPs in parallel for three scaffolds.

```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l mem=4G
#$ -j y
#$ -t 1-3
#$ -o gatk.log
#$ -N gatk

I=$(($SGE_TASK_ID-1))

# Regions (=scaffolds) that will be analysed
# There must be as many as tasks are specified above with '#$ -t'
REGIONS=(Hmel201001 Hmel201002 Hmel201003)

# Path to the genome
GENOME=/fastdata/$USER/varcal/genome/Hmel2.fa

# Path to directory with BAM files
ALIDIR=/fastdata/$USER/varcal/alignments_hisat2
# ALIDIR=/fastdata/$USER/varcal/alignments_rgs

# Path to output directory
OUTDIR=/fastdata/$USER/varcal/gatk

# Create output directory
mkdir -p $OUTDIR >& /dev/null

# Log for this region
LOG=$OUTDIR/gatk-${REGIONS[$I]}.log

# Redirect all the following screen outputs to log file
exec > $LOG 2>&1

hostname
date
echo "=============================================================================="

# Load genomics software repository
source /usr/local/extras/Genomics/.bashrc

# Create dictionary for genome
# gatk CreateSequenceDictionary -R $GENOME -O ${GENOME%.fa}.dict

cd $ALIDIR

# String for alignments
ALIGNMENTS=(*.bam)
BAMS=$(printf -- '-I %s ' "${ALIGNMENTS[@]}")

# Output vcf
VCF=$OUTDIR/gatk-${REGIONS[$I]}.vcf
# Output bcf
BCF=$OUTDIR/gatk-${REGIONS[$I]}.bcf

# Call SNPs and genotypes with GATK
# -----------------------------------------------------------
#  * --java-options "-Xmx4g": set the maximum memory for java to 4g 
#                             (should match what it is specified in
#                              the header #$ -l mem)
#   * HaplotypeCaller: use this caller
#   * --base-quality-score-threshold 20: filter out bases with mapping quality <20
#   * --minimum-mapping-quality 20: filter out alignments with mapping quality <20
#   * --output-mode EMIT_VARIANTS_ONLY: output only variants
#   * --dont-use-soft-clipped-bases: avoid using soft-clipped bases for calls

gatk --java-options "-Xmx4g" \
HaplotypeCaller \
--base-quality-score-threshold 20 \
--minimum-mapping-quality 20 \
--output-mode EMIT_VARIANTS_ONLY \
--dont-use-soft-clipped-bases \
-R $GENOME \
$BAMS \
-L ${REGIONS[$I]} \
-O $VCF >& $LOG

# convert from VCF to BCF
bcftools view -O b $VCF > $BCF

# index BCF file
bcftools index $BCF
# -----------------------------------------------------------

echo "=============================================================================="
date
```
We are going to use the tool called [``HaplotypeCaller``](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php). This is a more recent and improved caller with respect to the previous widely used ``UnifiedGenotyper``. The main improvements are a better identification of indels and the implementation of a multialleic model for >2 allele. We are using a number of options:
* `--java-options "-Xmx4g"`: GATK is a java program and as such it requires specifying the maximum memory that can be used
* `HaplotypeCaller`: use this caller
* `--base-quality-score-threshold 20`: filter out sites with base quality (BQ) <20
* `--minimum-mapping-quality 20`: filter out alignments with mapping quality (MQ) <20
* `--output-mode EMIT_VARIANTS_ONLY`: output only variants
* `--dont-use-soft-clipped-bases`: avoid using soft-clipped bases to minimize false positive and false negative calls; this is recommended for RNASeq (and it is how bcftools works).

When you have finished editing the bash script, save it as `gatk.sh`, make it executable with `chmod` and submit it to the job queue with `qsub`:
```bash
chmod +x gatk.sh
qsub gatk.sh
```
If all goes well, it should take no longer than a few minutes to get the jobs finished (although you will notice it takes a bit longer than bcftools). You should now have the following files:
```ls -lh gatk```

>``total 22M``<br>
>``-rw-r--r-- 1 myuser cs  92K Feb 16 14:28 gatk-Hmel201001.bcf``<br>
>``-rw-r--r-- 1 myuser cs  178 Feb 16 14:28 gatk-Hmel201001.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs 828K Feb 16 14:28 gatk-Hmel201001.log``<br>
>``-rw-r--r-- 1 myuser cs 334K Feb 16 14:28 gatk-Hmel201001.vcf``<br>
>``-rw-r--r-- 1 myuser cs  18K Feb 16 14:28 gatk-Hmel201001.vcf.idx``<br>
>``-rw-r--r-- 1 myuser cs  16M Feb 16 15:23 gatk-Hmel201002.log``<br>
>``-rw-r--r-- 1 myuser cs 5.0M Feb 16 15:21 gatk-Hmel201002.vcf``<br>
>``-rw-r--r-- 1 myuser cs  23K Feb 16 14:24 gatk-Hmel201003.bcf``<br>
>``-rw-r--r-- 1 myuser cs  157 Feb 16 14:24 gatk-Hmel201003.bcf.csi``<br>
>``-rw-r--r-- 1 myuser cs 378K Feb 16 14:24 gatk-Hmel201003.log``<br>
>``-rw-r--r-- 1 myuser cs  86K Feb 16 14:24 gatk-Hmel201003.vcf``<br>
>``-rw-r--r-- 1 myuser cs  18K Feb 16 14:24 gatk-Hmel201003.vcf.idx``<br>

The content of the logfiles can be quite long with many often harmless warnings. In our case, most of the warnings are due to missing data not allowing to calculate some VCF annotations. This can be annoying and result in pretty big files, but this can be avoiding by restricting the logging level to errors only with ``--verbosity ERROR``.

## 4. Operations with BCF files
The next sections exemplify how to do operations with VCF/BCF files, including merging, subsetting and filtering, mostly using bcftools and awk.

### Samples and SNPs 
A list of the samples contained in the file can be obtained using simple linux commands or `bcftools query`, and can be counted with `wc`:
```bash
bcftools view -h bcftools/bcftools-Hmel201001.bcf | grep "^#CHROM" | sed 's/.*FORMAT\t//; s/\t/\n/g'
bcftools view -h bcftools/bcftools-Hmel201001.bcf | grep "^#CHROM" | sed 's/.*FORMAT\t//; s/\t/\n/g' | wc -l
bcftools query -l bcftools/bcftools-Hmel201001.bcf
bcftools query -l bcftools/bcftools-Hmel201001.bcf | wc -l
```
The number of variants can be counted using `bcftools view -H` and `wc`, but it includes all kinds of variants. To check particular ones, we need to use the flags `-v/-V` (or `--types/--exclude-types`) to include or exclude certain variant. 
```bash
# all variants
bcftools view -H bcftools/bcftools-Hmel201001.bcf | wc -l
# only SNPs
bcftools view -H -v snps bcftools/bcftools-Hmel201001.bcf | wc -l
# only indels
bcftools view -H -v indel bcftools/bcftools-Hmel201001.bcf | wc -l
```
Multiallelic SNPs (with >2 alternate alleles) can be extracted with awk or with the `-m` flag:
```bash
bcftools view -H -v snps bcftools/bcftools-Hmel201001.bcf | awk '$5 ~ /,/' | less -S
bcftools view -H -v snps -m 3 bcftools/bcftools-Hmel201001.bcf | less -S
```

### Merging and subsetting files
Non-overlapping files with the same samples can be easily merged with `bcftools concat`. Let's merge the bcf files for our three scaffolds, index it and have a look:
```bash
bcftools concat bcftools/*.bcf -O b -o bcftools/bcftools-concat.bcf
bcftools index bcftools/bcftools-concat.bcf
bcftools view -H bcftools/bcftools-concat.bcf | less -S
```
Likewise, we can subset SNPs within particular regions, for example to extract the variants within the region 5000-5500bp in scaffold Hmel201001 and 15000-20000 in Hmel201002:
```bash
bcftools view -H bcftools/bcftools-concat.bcf Hmel201001:5000-5500,Hmel201002:15000-20000 | less -S 
```
Subsetting a number of samples is also possible, for example to get only the first 10 samples and the last 10 samples:
```bash
SAMPLES=$(bcftools query -l bcftools/bcftools-concat.bcf | head -n 10 | tr '\n/' ','| sed 's/,$//')
bcftools view -s $SAMPLES bcftools/bcftools-concat.bcf -O b -o bcftools/bcftools-concat-first16.bcf
bcftools index bcftools/bcftools-concat-first16.bcf

SAMPLES=$(bcftools query -lDepthPerSampleHC bcftools/bcftools-concat.bcf | tail -n 10 | tr '\n/' ','| sed 's/,$//')
bcftools view -s $SAMPLES bcftools/bcftools-concat.bcf -O b -o bcftools/bcftools-concat-last16.bcf
bcftools index bcftools/bcftools-concat-last16.bcf
```
And then merge them into a single file with `bcftools merge`:
```bash
bcftools merge bcftools/bcftools-concat-first16.bcf bcftools/bcftools-concat-last16.bcf -O b -o bcftools/bcftools-merge.bcf
bcftools index bcftools/bcftools-merge.bcf
```

### Comparing and combining outputs
We are now going to compare the SNPs called by bcftools and GATK. First, let's merge the bcf files for each scaffold:
``bash
mkdir comparison
# merge bcftools calls
bcftools concat bcftools/*.bcf -O b -o comparison/bcftools.bcf
bcftools index comparison/bcftools.bcf
# merge GATK calls
bcftools concat gatk/*.bcf -O b -o comparison/gatk.bcf
bcftools index comparison/gatk.bcf
```
Let's see the total number of SNPs called by each programme:
```bash
bcftools view -H comparison/bcftools.bcf | wc -l
bcftools view -H comparison/gatk.bcf | wc -l
```
As you can see both programmes called a similar number of SNPs. Let's now find out the SNPs that are shared between both calls using ``bcfools isec``:
```bash
bcftools isec -O b comparison/bcftools.bcf comparison/gatk.bcf -p comparison/isec
```
This will result in the following files:
```ls -lh comparison/isec```

>``total 1.7M``<br>
>``-rw-r--r-- 1 bo1vsx bo 1.5M Feb 16 16:20 0000.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo 1.6K Feb 16 16:20 0000.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo  35K Feb 16 16:20 0001.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo  217 Feb 16 16:20 0001.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo  70K Feb 16 16:20 0002.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo  229 Feb 16 16:20 0002.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo  80K Feb 16 16:20 0003.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo  226 Feb 16 16:20 0003.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo  559 Feb 16 16:20 README.txt``<br>

The content of the BCF files is explained in the README.txt file:
* ``0000.bcf``: private to bcftools
* ``0001.bcf``: private to gatk
* ``0002.bcf``: bcftools records shared with gatk
* ``0003.bcf``: gatk records shared with gatk

### Filtering
The steps above got us a quite raw sets of SNPs, but usually we want to apply some filters depending on the sort of downstream analyses that we intend to do. Filtering usually consist on establishin some sort of threshold for the SNP quality in the QUAL field and some of the variables encoded in the INFO field (e.g. depth, mapping quality, etc.), but more sophisticated filterin is possible based on genotype-specific variables (e.g. discarding genotypes below a certain quality or with fewer than a number of reads). Although in principle it is possible to do soft-filtering of SNPs and simply tag them as `PASSED` or `FAILED` in the VCF/BCF file, in practice it is more useful to do hard-filtering and simply remove all those SNPs (and genotypes) that did not meet our thresholds. We are going to focus on only the SNPs, so let's extract only the snps from the variants shared by bcftools and gatk:

```bash
mkdir filtering
bcftools view -v snps comparison/isec/0002.bcf -O b > filtering/snps.bcf
bcftools index filtering/snps.bcf
```

#### SNP-based filtering
bcftools allows applying filter either with ``bcftools view`` or with ``bcftools filter`` and using information encoded in the QUAL or INFO fields. These are some examples:

##### Filter by SNP quality - exclude SNPs with QUAL<20
```bash
bcftools view -e 'QUAL<20' -O b filtering/snps.bcf > filtering/snps.QS20.bcf
```

##### Filter by SNP depth - exclude SNPs with depth <30
```bash
bcftools view -e 'INFO/DP<100' -O b filtering/snps.bcf > filtering/snps.DP100.bcf
```
An important note here is to be aware that if the file is further processed so that only part of the individuals are used, the fields in INFO may not be updated and it would be then unreliable to filter out using any information from that field.

##### Filter by mapping quality

##### Filter by allele frequency
For population genetic analyses it is frequent to remove variants below a certain allele frequency, as these ones are difficult to tell apart from sequencing errors. For example, to exclude all SNPs with a minor allele frequency (MAF) below 5% we would run:
```bash
bcftools view -e 'MAF<0.05' -O b filtering/snps.bcf > filtering/snps.MAF005.bcf
```
##### Filter by sample coverage
Another typical situation is to want to exclude all SNPs for which only a small fracion of all the individuals have sequence data

#### Combining multiple filters

#### Genotype-based filtering
```bash
bcftools view -e ' SUM(FMT/DP)<100' -O b filtering/snps.bcf > filtering/snps.DP100.bcf
```
---
---
## 5. Extras

All the following steps are optional extras

### Population structure with NGSADMIX

### PCA of genoypes with R
We are going to carry out PCA using R. For that, we will need to convert our BCF file to mean genotype format (based on the BIMBAM format), where genotypes are encoded as mean genotypes. A mean genotype is a value between 0 and 2 that can be interpreted as the minor allele dosage: 0 is homozygous for the major allele, 1 is a heterozygote, and 2 is a homozygote for the minor allele. Thus, intermediate values reflect the uncertainty in genotypes. 

We are going to use a custom Perl script called ``bcf2bbgeno.pl`` to get such mean genotypes. This scripts (1) removes all the SNPs that have more than two alleles and (2) calculates empirical posterior genotype probabilities from the genotype likelihoods in the BCF file under the assumption that the population is in Hardy-Weinberg equilibrium (HWE). Specifically, the script uses inferred allele frequencies to set HWE priors:

p(AA) = p2; p(aa) = (1-p)2; p(Aa) = 2p(1-p)

being p the allele frequency of major/reference allele A. Genotype likelihoods are multiplied by these priors to obtain genotype posterior probabilities that are then encoded as mean genotypes and saved to a .bbgeno file.

You can get some info about how to run the Perl script:
```bash
mkdir scripts
# download script
wget https://raw.githubusercontent.com/visoca/popgenomworkshop-gwas_gemma/master/scripts/bcf2bbgeno.pl -O scripts/bcf2bbgeno.pl
# give execution permissions
chmod +x scripts/bcf2bbgeno.pl
# show help
scripts/bcf2bbgeno.pl -h
```

To calculate the genotype posterior probabilites and save them in mean genotype format, we would need to run a command like this one, including the flag to include a header with the individuals names (followed by compression to save some space, given that ```R``` can handle gzipped files without problems):
```bash
scripts/bcf2bbgeno.pl -i snps/snps.flt.bcf -o snps.bbgeno -p H-W -a
# compress the file to save some space
gzip fha.bbgeno
```
We will also need the information about the race and sex of the samples, which can be copied from a shared directory in ShARC:
```bash
cp
```
Then we will use R to do PCA with the mean genotypes and do some plotting. First, let's load the genotypes
```R
## load genotypes
genofile<-"pca/snps.bbgeno.gz"
genotypes<-read.table(genofile,header=T, check.names=F)
```
Now we will calculate the genotype covariance matrix:
```R
## calculate N x N genotype covariance matrix
gmn<-apply(genotypes[,-(1:3)],1,mean)
nids<-ncol(genotypes)-3
nloci<-nrow(genotypes)
gmnmat<-matrix(gmn,nrow=nloci,ncol=nids)
gprime<-genotypes[,-(1:3)]-gmnmat

gcovarmat<-matrix(NA,nrow=nids,ncol=nids)
for(i in 1:nids){
 for(j in i:nids){
  if (i==j){
    gcovarmat[i,j]<-cov(gprime[,i],gprime[,j])
  }
  else{
    gcovarmat[i,j]<-cov(gprime[,i],gprime[,j])
    gcovarmat[j,i]<-gcovarmat[i,j]
  }
 }
}
```
And we run a PCA using the function `prcomp`:
```R
## pca on the genotype covariance matrix
pcg<-prcomp(x=gcovarmat,center=TRUE,scale=FALSE)
```
Extract the PCS:
```R
pcs<-pcg$x
rownames(pcs)<-colnames(genotypes[,-(1:3)])
pcs.col<-unlist(lapply (rownames(pcs), function(x)
as.character(col.table[col.table[,2]==populations[populations[,1]==as.character(x),][[2]],][[1]])))
```
```R

# plot each PC using different colours for race and different symbols for sex
race.col<-c("","")
sex.symbol<-c(0,2)
#Plot PC1 vs PC2
plot(pcs[,1],pcs[, i+2], main="PC1 vs PC2", xlab="PC1", ylab="PC2"))
```

We can also plot the names of the samples to investigate potential outliers:
```R
# plot each PC with sample names
for (i in 1:(npcs-1)){
par(mfrow=c(1,1),oma=c(1,1,1,10), mar=c(5, 4, 4, 2) + 0.1)
plot(pcs[,i],pcs[,i+1],type="n", main="PCA using covariance matrix", xlab=paste("PC",i,sep=""), ylab=paste("PC",(i+1),sep=""))
text(pcs[,i],pcs[,i+1],labels=rownames(pcs),col=pcs.col,cex=0.5)
# plot legend to the right
par(mfrow=c(1,1), oma=c(0,0,0,1), mar=c(0,0,0,0), new=T)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
legend("right", legend=pops, lty=1, lwd=10, col=colours, bty="n", xpd=T)

```
### Differentiation genome scans using pairwise FSTs

### SNP and genotype calling with ANGSD

### Resources
* [samtools manual](http://www.htslib.org/doc/samtools.html)
* [bcftools manual](http://www.htslib.org/doc/bcftools.html)
* [bcftools howto](http://samtools.github.io/bcftools/howtos/index.html)
* [GATK User Guide](https://software.broadinstitute.org/gatk/documentation/quickstart?v=4)
* [ANGSD manual](http://www.popgen.dk/angsd/index.php/ANGSD)
* [SAM/BAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
* [VCF/BCF format specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
* [VCF/BCF format visual explanation](http://vcftools.sourceforge.net/VCF-poster.pdf) - __Highly recommended__
 
### References
* [Li 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/) - bcftools mpileup implementation
* [Danacek et al 2014](http://samtools.github.io/bcftools/call-m.pdf) - bcftools multiallelic caller
* [Van der Auwera et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/) - GATK Best Practices for Variant Discovery for beginners. Note some procedures may be out-dated, check the current documentation.
* [Novembre et al 2008](https://www.nature.com/articles/nature07331) - Example of using PCAs of genotypes to investigate population structure
* [Korneliussen et al. 2014](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4) - ANGSD publication
* [Bhatia et al. 2013](http://genome.cshlp.org/content/23/9/1514.full) - Excellent paper about F<sub>ST</sub> estimation and interpretation.

Add key Heliconius references??
