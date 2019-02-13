**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# SNP and genotype calling
#### Victor Soria-Carrasco

The aim of this practical is to learn how to call single nucleotide polymorphism (SNPs) and genotypes, that is the process of identifying variable sites and determining the genotype for each individual at each site. We will be using a dataset of whole genome sequence data of 32 individuals of *Heliconius melpomene*. After calling SNPs, we will do some subsetting and filtering and will carry out a few example analyses 

### Resources
* [samtools manual](http://www.htslib.org/doc/samtools.html)
* [bcftools manual](http://www.htslib.org/doc/bcftools.html)
* [bcftools howto](http://samtools.github.io/bcftools/howtos/index.html)
* [GATK User Guide](https://software.broadinstitute.org/gatk/documentation/quickstart?v=4)
* [Van der Auwera et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/) - GATK Best Practices for Variant Discovery for beginners. Note some procedures may be out-dated, check the current documentation.
* [ANGSD manual](http://www.popgen.dk/angsd/index.php/ANGSD)
* PCA ref
* [Korneliussen et al. 2014](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4) ANGSD publication
* [Bhatia et al. 2013](http://genome.cshlp.org/content/23/9/1514.full) - Excellent paper about F<sub>ST</sub> estimation and interpretation.
Add key Heliconius references??

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
We are going to use three different tools for SNP calling: [bcftools](http://www.htslib.org/), [GATK](https://software.broadinstitute.org/gatk/), and, optionally,[ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). Additionally, we will use bcftools and awk for filtering files. All of them are already installed as part of the [Genomics Software Repository](http://soria-carrasco.staff.shef.ac.uk/softrepo/).

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

Lastly, we also need to download the reference genome:
```bash
mkdir genome
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2_-_scaffolds.fa.gz -O genome/Hmel2.fa.gz
```

## BCFtools SNP calling
BCFtools is a very popular programme to call SNPs and genotypes (and also to manipulate and filter vcf/bcf files as we will see afterwards). SNP calling is a relatively intensive process, to speed things up we will be restricting variant calling to 3 scaffolds. Before calling SNPs, we have to index the genome using samtools. We can keep the genome compressed to save space - which is always a good idea when working with NGS data -, but we have to use bgzip (part of samtools) instead of gzip:
```bash
gzip -d genome/Hmel2.fa.gz
bgzip genome/Hmel2.fa
samtools faidx genome/Hmel2.fa.gz
```
This will produce two index files:
```bash
 ls -lh genome
```
>``-rw-r--r-- 1 myuser cs 76M Feb 13 17:11 Hmel2.fa.gz``<br>
>``-rw-r--r-- 1 myuser cs 31K Feb 13 17:11 Hmel2.fa.gz.fai``<br>
>``-rw-r--r-- 1 myuser cs 66K Feb 13 17:11 Hmel2.fa.gz.gzi``<br>

We are going to prepare now a batch script to call SNPs and genotypes for all 32 individuals. To speed things up, we will be using an [SGE array job](http://docs.hpc.shef.ac.uk/en/latest/parallel/JobArray.html) to call SNPs for three scaffolds (Hmel201001, Hmel201002, and Hmel201003) in parallel:

```bash {.line-numbers}
#!/bin/bash
#$ -l h_rt=2:00:00
#$ -l mem=2G
#$ -t 1-3
#$ -m bea
#$ -M myemail@mail.com
#$ -j y
#$ -o bcftools.log
#$ -N bcftools

# internal SGE task index
I=$(($SGE_TASK_ID-1))

hostname
date
echo "=============================================================================="

# Load genomics software repository
# (safer to load it here in case the worker nodes don't inherit the environment)
source /usr/local/extras/Genomics/.bashrc 

# Path to the genome
GENOME=/fastdata/$USER/varcal/genome/Hmel2.fa.gz
# Path to directory with BAM files
ALIDIR=/fastdata/$USER/varcal/alignments
# Path to output directory
OUTDIR=/fastdata/$USER/varcal/bcftools

# Regions (=scaffolds) that will be analysed
# There must be as many as tasks are specified above with '#$ -t'
REGIONS=(Hmel201001 Hmel201002 Hmel201003)

# Create output directory
mkdir -p $OUTDIR 

# Change to alignment dir
cd $ALIDIR

# Call SNPs and genotypes with bcftools
#  mpileup
#   * -q 20: filter out alignments with mapping quality <20
#   * -Q 20: filter bases with QS < 20
#   * -P ILLUMINA: restrict to Illumina platform for indels
#  call
#   * -m: use the multiallelic caller
#   * -v: output variants only
#   * -P 1e-6: prior on mutation rate

bcftools mpileup -Ou \
--max-depth 10000 \
-q 20 -Q 20 -P ILLUMINA \
-f $GENOME \
-r ${REGIONS[$I]} *.bam | \
bcftools call -mv \
-P 1e-6 \
-O b \
-o $OUTDIR/bcftools-${REGIONS[$I]}.bcf

# Index bcf file
bcftools index $OUTDIR/bcftools-${REGIONS[$I]}.bcf

echo "=============================================================================="
date
```
Calling SNPs with bcftools is a two-step process. First, `bcftools mpileup` estimate genotype likelihoods at each genomic position with sequence data. Second, `bcftools call` do identify both variants and genotypes, i.e. makes the actual call. To avoid generating intermediate temporary files, the output of `bcftools mpileup` is *piped* to `bcftools call`. We are also filtering out sites with base quality (BQ) <20 (-q 20) and alignemnts with mapping quality (MQ) <20 (-Q 20) and restrict the indel to ILLUMINA platform. We will use the multiallelic caller (-m), which allow more than two alleles, output only variant (-v), and use a relatively stringent prior for the expected subsitution rate (-P 1e-6). Lastly, we are saving the output in BCF format, which is the binary version of VCF (more details below).

## VCF and BCF format
The most popular text-file format for storing genetic variation information is the [Variant Call Format (VCF)](http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)). There are different versions of the format, but the core elements are the same. You can find a full description of the latest iteration [here](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf).Due to the large amount of data usually involved, files tend to be stored compressed. BCF is the binary version of VCF, which keeps the same information, but compressed and indexed. It  is designed to be much more efficient to process and store large amounts of data. The practical downside is that the contents can only be accessed with `bcftools view`.

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
10 to end | Sample id (from SM tag in BAM file RG string) | Sample genotype information (Usually also stores the genotype qualities (GQ) and genotype likelihoods (GL or PL (phred-scaled)
 
Let's now have a look at our files:
```bash
bcftools view varcal.bcf
```

## Operations with VCF/BCF files

When the job has finished, we will proceed to merge all the bcf files:


## GATK SNP calling

## ANGSD SNP calling (optional)


## Subsetting SNPs
### Extracting a set of SNPs for a region
### Extracting a set of SNPs for a set of individuals

## Filtering SNPs

## Comparing outputs

## Extras
### Population structure with NGSADMIX
### PCA of genoypes with R
### Quick genome scan using pairwise FSTs

