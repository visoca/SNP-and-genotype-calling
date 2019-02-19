**_Victor Soria-Carrasco_**

**_Department of Animal and Plant Sciences, University of Sheffield_**

This practical is part of the module [Advanced Data Analysis - Introduction to NGS data analysis](https://visoca.github.io/AdvDataAna-introNGS/). The aim is to learn how to call single nucleotide polymorphism (SNPs) and genotypes, that is the process of identifying variable sites and determining the genotype for each individual at each site. We will be using a dataset of whole genome sequence data of 32 individuals of *Heliconius melpomene*. After calling SNPs, we will do some subsetting and filtering and will carry out a few example analyses.

## Table of contents
1. [Initial set up](setup.md)
2. [SNP and genotype calling with BCFtools](bcftools.md)
3. [VCF and BCF format](vcfbcf.md)
4. [SNP and genotype calling with GATK](gatk.md)
5. [Operations with BCF files](bcfops.md)

Extras:<br>
  * [Population structure with NGSADMIX](ngsadmix.md)<br>
  * [PCA of genoypes with R](pca.md)<br>
  * [SNP and genotype calling with ANGSD](angsd.md)<br>
  * [Genetic architecture of traits with GEMMA](https://visoca.github.io/popgenomworkshop-gwas_gemma/)<br>
  * [Delimitation of contiguous regions of differentiation using Hidden Markov Models](https://visoca.github.io/popgenomworkshop-hmm/)<br>

---

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
