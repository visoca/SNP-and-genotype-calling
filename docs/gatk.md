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

Now, let's prepare a batch script to run GATK on ShARC. As before, we will use an SGE array to call SNPs in parallel for three scaffolds, using 4 cores for each task.

```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l mem=4G
#$ -j y
#$ -t 1-3
#$ -pe smp 4
#$ -o gatk.log
#$ -N gatk

I=$(($SGE_TASK_ID-1))

# Regions (=scaffolds) that will be analysed
# There must be as many as tasks are specified above with '#$ -t'
REGIONS=(Hmel201001 Hmel201003 Hmel201008)

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
#   * --native-pair-hmm-threads 4: use 4 threads for the pair-HMM process

gatk --java-options "-Xmx4g" \
HaplotypeCaller \
--base-quality-score-threshold 20 \
--minimum-mapping-quality 20 \
--output-mode EMIT_VARIANTS_ONLY \
--dont-use-soft-clipped-bases \
--native-pair-hmm-threads 4 \
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
* --native-pair-hmm-threads 4: use 4 threads to speed up the pair-HMM process

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

[Back to TOC](index.md)
