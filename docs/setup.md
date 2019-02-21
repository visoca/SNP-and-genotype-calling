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
We need to copy alignments in BAM format. We are not going to use the same alignments used in the previous sessions, because we need them to contain certain specific information to run some of the exercises for this practical. You can copy them with rsync, which allows resuming in case something goes wrong and the transfer is interrupted:
```bash
rsync -avPh /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/SNPgenocall/alignments ./
```
If the transfer is too slow, you can cancel it (ctrl+c) and try from a different location in `fastdata`:
```bash
rsync -avPh /fastdata/SNPgenocall/alignments ./
```
It is important to have all the BAM files indexed. To speed things up, we can use an [SGE array job](http://docs.hpc.shef.ac.uk/en/latest/parallel/JobArray.html) to index files in parallel:
```bash
#!/bin/bash
#$ -l mem=2g
#$ -l h_rt=1:00:00
#$ -t 1-32
#$ -j y
#$ -o indexbam.log
#$ -N indexbam

# Load genomics software repository
# (safer to load it here in case the worker nodes don't inherit the environment)
source /usr/local/extras/Genomics/.bashrc

I=$(($SGE_TASK_ID-1))

BAMS=($(ls alignments/*.bam))

samtools index ${BAMS[$I]}
```
When you have finished editing the bash script, save it as `indexbam.sh`, make it executable with `chmod` and submit it to the job queue with `qsub`:
```bash
chmod +x indexbam.sh
qsub indexbam.sh
```
Indexing should be finished in a few minutes and the files should look like this:
```bash
ls -lh alignments | head -n 8
```
>``total 59G``<br>
>``-rwx------ 1 myuser cs 1.9G Feb 14 14:46 60A.bam``<br>
>``-rw-r--r-- 1 myuser cs 486K Feb 15 03:18 60A.bam.bai``<br>
>``-rwx------ 1 myuser cs 1.8G Feb 14 14:44 60I.bam``<br>
>``-rw-r--r-- 1 myuser cs 471K Feb 15 03:18 60I.bam.bai``<br>
>``-rwx------ 1 myuser cs 2.1G Feb 14 14:39 61A.bam``<br>
>``-rw-r--r-- 1 myuser cs 590K Feb 15 03:18 61A.bam.bai``<br>
>``-rwx------ 1 myuser cs 1.9G Feb 14 14:41 61I.bam``<br>
>``-rw-r--r-- 1 myuser cs 505K Feb 15 03:18 61I.bam.bai``<br>

We also need to download the reference genome, you can do it using [wget](https://www.gnu.org/software/wget/manual/wget.html) (remember that you can also read the manual on the command line with `man wget`):
```bash
mkdir genome
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2_-_scaffolds.fa.gz -O genome/Hmel2.fa.gz
```
[Back to TOC](index.md)
