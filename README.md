# SNP and genotype calling

We will be using a dataset of X individuals of Heliconius melpomene for this SNP calling tutorial. We will learn how to identify variants and 

## Initial set up
We will be using the University of Sheffield HPC cluster [ShARC](https://www.sheffield.ac.uk/cics/research/hpc/sharc) and several programmes installed as part of the [Genomics Software Repository](http://soria-carrasco.staff.shef.ac.uk/softrepo/). Please ensure you have set up your account to use the repository. If the repository is set up correctly, you should see the following message every time you call an interactive job in a working node:
```
  Your account is set up to use the Genomics Software Repository
    More info: http://soria-carrasco.staff.shef.ac.uk/softrepo
```
If you are interested in knowing what software is currently installed, you can use the command `softrepo`.

For this particular tutorial, it is expected that you will be working on a directory called `varcal` in your /fastdata/$USER directory:
```bash
cd /fastdata/$USER
mkdir varcal
cd varcal
```
Remember you can check your current working directory anytime with the command `pwd`.
It should show something like:
```bash
pwd
/fastdata/bo11zzz/varcal
```

## Programmes
We are going to use three different tools for SNP calling: [samtools/bcftools](http://www.htslib.org/), [GATK](https://software.broadinstitute.org/gatk/), and [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). Additionally, we will use bcftools and awk for filtering files. All of them are already installed as part of the [Genomics Software Repository](http://soria-carrasco.staff.shef.ac.uk/softrepo/).

## Data

## samtools/bcftools SNP calling
To speed things up, we will be restricting variant calling to 10 scaffolds that cover X bp.
## GATK SNP calling

## ANGSD SNP calling

## VCF and BCF formats
The most popular text-file format for storing genetic variation information is the [Variant Call Format (VCF)](http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)). There are different versions of the format, but the core elements are the same. You can find a full description of the latest iteration [here](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf).Due to the large amount of data usually involved, files tend to be stored compressed. 

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

## Subsetting SNPs
### Extracting a set of SNPs for a region
### Extracting a set of SNPs for a set of individuals

## Filtering SNPs

## Comparing outputs

## Extras
### Population structure with NGSADMIX
### PCA of genoypes with R
### Quick genome scan using pairwise FSTs

