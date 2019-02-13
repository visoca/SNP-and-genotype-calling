**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# SNP and genotype calling
#### Victor Soria-Carrasco

The aim of this practical is to learn how to call single nucleotide polymorphism (SNPs) and genotypes, that is the process of identifying variable sites and determining the genotype for each individual at each site. We will be using a dataset of whole genome sequence data of 32 individuals of *Heliconius melpomene*. After calling SNPs, we will do some subsetting and filtering and will carry out a few example analyses 

---
layout: "Resources"
thumb: string
title: string
---
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
/fastdata/bo11xxx/varcal
```

## Programmes
We are going to use three different tools for SNP calling: [bcftools](http://www.htslib.org/), [GATK](https://software.broadinstitute.org/gatk/), and, optionally,[ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). Additionally, we will use bcftools and awk for filtering files. All of them are already installed as part of the [Genomics Software Repository](http://soria-carrasco.staff.shef.ac.uk/softrepo/).

## Data
We will need to copy the reference genome and the alignments in BAM format produced in the previous sessions to `/fastdata/$USER/varcal/alignments`. If you don't have such data, there are files are available in ` /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/data`, e.g.:
```bash
cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/data/genome /fastdata/$USER/varcal/
cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/data/alignments /fastdata/$USER/varcal/
```


## BCFtools SNP calling
BCFtools is a very popular programme to call SNPs and genotypes (and also to manipulate and filter vcf/bcf files as we will see afterwards). SNP calling is a relatively intensive process, to speed things up we will be restricting variant calling to 10 scaffolds that cover X bp. Before calling SNPs, we have to index the genome using
When the job has finished, we will proceed to merge all the bcf files:
```bash
  
```

## GATK SNP calling

## ANGSD SNP calling (optional)

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

