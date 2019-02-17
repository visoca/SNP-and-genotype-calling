## Population structure with NGSadmix

[NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) is a method for inferring admixture from genotype likelihoods, thus taking into account the uncertainty inherent to NGS data. NGSadmix is part of the package ANGSD, it uses genotype likelihood information rather than absolute genotype calls. 

We first need to conver the bcf file to beagle format with angsd:
```bash
mkdir ngsadmix
angsd -vcf-gl <(bcftools view -O v filtering/snps.bcf) -fai genome/Hmel2.fa.fai -doMaf 3 -nInd 32 -domajorminor 1 -doglf 2 -out ngsadmix/snps
```


