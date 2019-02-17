
## 2. VCF and BCF format
One of the most popular text-file formats for storing genetic variation information is the [Variant Call Format (VCF)](http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)). There are different versions of the format, but the core elements are the same. You can find a full description of the latest iteration [here](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf).Due to the large amount of data usually involved, files tend to be stored compressed. BCF is the binary version of VCF, which keeps the same information, but compressed and indexed. It  is designed to be much more efficient to process and store large amounts of data. The practical downside is that the contents can only be accessed with `bcftools view`.

The VCF format is composed of meta-information lines (prefixed wih ##), a header line (prefixed with #), and data lines each containing information about a position in the genome and genotype information on samples for each position (text fields separated by tabs). There are 8 fixed field lines, usually followed by an additional FORMAT field and an arbitrary number of sample fields when genotypes are included.

**Table 1** Description of VCF fields

Column Number| Title | Description
:--|:--|:--
1 | CHROM | Chromosome/scaffold of the reference genome
2 | POS | Position on the chromosome/scaffold given in CHROM (1-based)
3 | ID | ID for the SNP/INDEL (Here always '.')
4 | REF | Reference allele (usually the allele on reference genome). I can have several bases in the case of indels.
5 | ALT | Alternate alleles (separated by commas), can have several bases in the case of indels
6 | QUAL | Variant quality score ([phred-scale](https://gatkforums.broadinstitute.org/gatk/discussion/4260/how-should-i-interpret-phred-scaled-quality-scores))
7 | FILTER | Filter status ('.' if unfiltered or hard-filtered, PASS or a semicolon-separated list of codes for filters that fail if soft-filtered)
8 | INFO | List of variant annotations delimited by ';' and as specified in the meta-information lines
9 | FORMAT | Describes the format for the sample genotype information in column 10 onward
10 to end | Sample id (from SM tag in BAM file RG string) | Sample genotype information. Usually it stores the genotype call (GT), genotype likelihoods (GL or PL (phred-scaled)), depth (DP), allelic depth (AD), and genotype quality (GQ).

Some notes about genotype fields:
* `GT`: Allele values are 0 for the reference and 1+ alternate alleles. Allele are separated with `/` if unphased and with `|` if phased. There are as many alleles separated by `/`  or `|` as the ploidy. A `.` is used for missing data, i.e. when a genotype could not be called.
* `GL` or `PL`: Genotype likelihoods as floating-point log10-scaled or phred scaled format. A diploid biallelic has three values: homozygote for reference, heterozygote, and homozygote for alternate. Triallelic sites have 6 values: homozygote for reference, heterozygote of reference and alternate allele 1, homozygote for alternate allele 1, heterozygote of reference and alternate allele 2, heterozygote of alternate alleles 1 and 2, and homozygote for alternate allele 2.

A graphical explanation of the different parts and fields can be found [here](http://vcftools.sourceforge.net/VCF-poster.pdf).

Let's now have a look at our files:
```bash
bcftools view bcftools/Hmel201001.bcf | less -S
```
The file starts with a few lines specifying the VCF and bcftools versions, the command used with `bcftools mpileup`, followed by a long list of the >700 scaffolds that comprise the genome assembly (starting with `##contig`), the descriptions of the fields for ALT, INFO, and FORMAT, and any additional commands that were executed to produce this file (`bcftools call` in this case).

To focus only on the header or the SNPs, the flags `-h` and `-H` can be used:
```bash
bcftools view -h bcftools/Hmel201001.bcf | less -S
bcftools view -H bcftools/Hmel201001.bcf | less -S
```

[Back to TOC](index.md)
