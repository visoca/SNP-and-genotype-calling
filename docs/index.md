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
