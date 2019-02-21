## 1. Initial set up
First of all, please remind that this tutorial must be run using an interactive session in ShARC. For that, you should log in into ShARC with `ssh`, and then request an interactive session with `qrsh`. Your shell prompt should show `sharc-nodeXXX` (XXX being a number between 001 and 172) and not `@sharc-login1` nor `@sharc-login2`.

For this particular practical, we are going to create and work on a directory called `varcal` in your /fastdata/$USER directory:
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
There are many different tools for SNP calling, we are going to use three of the most popular: [bcftools](http://www.htslib.org/), [GATK](https://software.broadinstitute.org/gatk/), and, optionally, [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). Additionally, we will use bcftools and awk for manipulating VCF/BCF files. Additional exercises will use [R](https://cran.r-project.org/), [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) (part of ANGSD), and [gemma](http://www.xzlab.org/software.html). All of them are already installed as part of the [Genomics Software Repository](http://soria-carrasco.staff.shef.ac.uk/softrepo/).

## Data
We need to copy alignments in BAM format. We are not going to use the same alignments used in the previous sessions, because we need them to contain certain specific information to run some of the exercises for this practical. You can copy them with:
```bash
cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/SNPgenocall/alignments_hisat2 ./alignments
```
It is important that all the BAM files are indexed, which can be done with the following commands:
```bash
ls alignments/*.bam | xargs -I {} sh -c 'samtools index {}'
```
The files should look like this:
```bash
ls -lh alignments
```
>``-rw-r--r--   1 myuser cs 1.9G Feb 13 04:29 60A.bam``<br>
>``-rw-r--r--   1 myuser cs 482K Feb 13 04:30 60A.bam.bai``<br>
>``-rw-r--r--   1 myuser cs 1.7G Feb 13 04:29 60I.bam``<br>
>``-rw-r--r--   1 myuser cs 468K Feb 13 04:29 60I.bam.bai``<br>
>``-rw-r--r--   1 myuser cs 2.1G Feb 13 04:31 61A.bam``<br>
>``-rw-r--r--   1 myuser cs 586K Feb 13 04:31 61A.bam.bai``<br>

Lastly, we also need to download the reference genome, you can do it using [wget](https://www.gnu.org/software/wget/manual/wget.html) (remember that you can also read the manual on the command line with `man wget`):
```bash
mkdir genome
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2_-_scaffolds.fa.gz -O genome/Hmel2.fa.gz
```
[Back to TOC](index.md)
