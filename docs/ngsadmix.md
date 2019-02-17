## Population structure with NGSadmix

[NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) is a method for inferring admixture from genotype likelihoods, thus taking into account the uncertainty inherent to NGS data. NGSadmix is part of the package ANGSD, it uses genotype likelihood information rather than absolute genotype calls. 

We first need to conver the bcf file to beagle format with angsd:
```bash
mkdir ngsadmix
angsd -vcf-gl <(bcftools view -O v filtering/snps.bcf) -fai genome/Hmel2.fa.fai -doMaf 3 -nInd 32 -domajorminor 1 -doglf 2 -out ngsadmix/snps
```

Now prepare the script to run NGSadmix for K=2 to K=4 in parallel:

```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l mem=2G
#$ -t 1-4
#$ -m bea
#$ -M myemail@mail.com
#$ -j y
#$ -o ngsadmix.log
#$ -N ngsadmix

# internal SGE task index
K=$(($SGE_TASK_ID+1))


hostname
date
echo "=============================================================================="

NGSadmix -likes ngsadmix/snps.beagle.gz -K $K -o ngsadmix/snps_$K >& ngsadmix/snps_$K.sclog

echo "=============================================================================="
date
```
When you have finished editing the bash script, save it as ngsadmix.sh, make it executable with chmod and submit it to the job queue with qsub:
```bash
chmod +x ngsadmix.sh
qsub ngsadmix.sh
```
After running, you should have the following files: 
```bash
ls -lh ngsadmix
```

Now we are going to visualize the admixture estimates using ´R´.

```R
ls -lh ngsadmix
```

