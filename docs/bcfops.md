## 4. Operations with BCF files
The next sections exemplify how to do operations with VCF/BCF files, including merging, subsetting and filtering, mostly using bcftools and awk.

### Samples and SNPs 
A list of the samples contained in the file can be obtained using simple linux commands or `bcftools query`, and can be counted with `wc`:
```bash
bcftools view -h bcftools/bcftools-Hmel201001.bcf | grep "^#CHROM" | sed 's/.*FORMAT\t//; s/\t/\n/g'
bcftools view -h bcftools/bcftools-Hmel201001.bcf | grep "^#CHROM" | sed 's/.*FORMAT\t//; s/\t/\n/g' | wc -l
bcftools query -l bcftools/bcftools-Hmel201001.bcf
bcftools query -l bcftools/bcftools-Hmel201001.bcf | wc -l
```
The number of variants can be counted using `bcftools view -H` and `wc`, but it includes all kinds of variants. To check particular ones, we need to use the flags `-v/-V` (or `--types/--exclude-types`) to include or exclude certain variant. 
```bash
# all variants
bcftools view -H bcftools/bcftools-Hmel201001.bcf | wc -l
# only SNPs
bcftools view -H -v snps bcftools/bcftools-Hmel201001.bcf | wc -l
# only indels
bcftools view -H -v indels bcftools/bcftools-Hmel201001.bcf | wc -l
```
Multiallelic SNPs (with >2 alternate alleles) can be extracted with awk or with the `-m` flag:
```bash 
bcftools view -H -v snps bcftools/bcftools-Hmel201001.bcf | awk '$5 ~ /,/' | less -S
bcftools view -H -v snps -m 3 bcftools/bcftools-Hmel201001.bcf | less -S
```

### Extracting information
Partial information can be extracted using the `bcftools query`. In the example above we saw how to get the list of samples using the `l` option, but it can also be used to extract any fields using the `-f` option. For example, you can simple extract the list of positions with:
```bash 
bcftools query -f '%POS\n'  bcftools/bcftools-Hmel201001.bcf | less -S
```
Or combine multiple fields in the output, for example, scaffold/chromosome, position, reference and alternate alleles:
```bash 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n'  bcftools/bcftools-Hmel201001.bcf | less -S
```
You can also extract genotypes (notice the position of the tab to avoid adding an extra tab at the end of the lines):
```bash 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'  bcftools/bcftools-Hmel201001.bcf | less -S
```
Or genotype depths:
```bash 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP\t]\n'  bcftools/bcftools-Hmel201001.bcf | less -S
```
Or allele depths:
```bash 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n'  bcftools/bcftools-Hmel201001.bcf | less -S
```
Or genotype likelihoods:
```bash 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%PL]\n'  bcftools/bcftools-Hmel201001.bcf | less -S
```
You can even combine this with awk, Perl or other tools to do simple calculations. For example, you can calculate allele frequencies from the AC and AN counts:
```bash 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' bcftools/bcftools-Hmel201001.bcf | less -S
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' bcftools/bcftools-Hmel201001.bcf| awk '{print $1,$2,$3,$4,$5,$6,$5/$6}' | less -S
```
Or the mean genotype quality:
```bash 
# For all genotypes (uncalled ones = 0)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GQ]\n' bcftools/bcftools-Hmel201001.bcf | \
awk '{SUM=0; for(i=5; i<=NF; i++){SUM+=$i}; AVG=SUM/NF; print $1,$2,$3,$4,AVG}' | less -S
# Only for called genotypes:
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GQ]\n' bcftools/bcftools-Hmel201001.bcf | \
awk '{SUM=0; N=0; for(i=5; i<=NF; i++){SUM+=$i; if ($i>0) N+=1}; AVG=SUM/N; print $1,$2,$3,$4,AVG}' | less -S
```
### Merging and subsetting files
Non-overlapping files with the same samples can be easily merged with `bcftools concat`. Let's merge the bcf files for our three scaffolds, index it and have a look:
```bash
bcftools concat bcftools/*.bcf -O b -o bcftools/bcftools-concat.bcf
bcftools index bcftools/bcftools-concat.bcf
# Now we have more SNPs
bcftools view -H bcftools/bcftools-concat.bcf | wc -l
# And you can see other scaffolds as you scroll down
bcftools view -H bcftools/bcftools-concat.bcf | less -S
```
Likewise, we can subset SNPs within particular regions, for example, we can extract all the SNPs on a particular chromosome/scaffold:
```bash
bcftools view -H bcftools/bcftools-concat.bcf Hmel201008 | wc -l
bcftools view -H bcftools/bcftools-concat.bcf Hmel201008 | less -S 
```
Or extract the variants within the region 5000-5500bp in scaffold Hmel201001 and 15000-20000 in Hmel201003:
```bash
bcftools view -H bcftools/bcftools-concat.bcf Hmel201001:5000-5500,Hmel201003:15000-20000 | wc -l
bcftools view -H bcftools/bcftools-concat.bcf Hmel201001:5000-5500,Hmel201003:15000-20000 | less -S 
```
Subsetting a number of samples is also possible, for example to get only the first 10 samples and the last 10 samples:
```bash
SAMPLES=$(bcftools query -l bcftools/bcftools-concat.bcf | head -n 10 | tr '\n/' ','| sed 's/,$//')
bcftools view -s $SAMPLES bcftools/bcftools-concat.bcf -O b -o bcftools/bcftools-concat-first10.bcf
bcftools index bcftools/bcftools-concat-first10.bcf

SAMPLES=$(bcftools query -l bcftools/bcftools-concat.bcf | tail -n 10 | tr '\n/' ','| sed 's/,$//')
bcftools view -s $SAMPLES bcftools/bcftools-concat.bcf -O b -o bcftools/bcftools-concat-last10.bcf
bcftools index bcftools/bcftools-concat-last10.bcf
```
And then merge them into a single file with `bcftools merge`:
```bash
bcftools merge bcftools/bcftools-concat-first10.bcf bcftools/bcftools-concat-last10.bcf -O b -o bcftools/bcftools-merge20.bcf
bcftools index bcftools/bcftools-merge20.bcf
```

### Comparing and combining outputs
We are now going to compare the SNPs called by bcftools and GATK. First, let's merge the bcf files for each scaffold:
```bash
mkdir comparison
# merge bcftools calls
bcftools concat bcftools/bcftools-Hmel*.bcf -O b -o comparison/bcftools.bcf
bcftools index comparison/bcftools.bcf
# merge GATK calls
bcftools concat gatk/gatk-Hmel*.bcf -O b -o comparison/gatk.bcf
bcftools index comparison/gatk.bcf
```
Let's see the total number of SNPs called by each programme:
```bash
bcftools view -H comparison/bcftools.bcf | wc -l
bcftools view -H comparison/gatk.bcf | wc -l
```
As you can see both programmes called a similar number of SNPs. Let's now find out the SNPs that are shared between both calls using ``bcfools isec``:
```bash
bcftools isec -O b comparison/bcftools.bcf comparison/gatk.bcf -p comparison/isec
```
This will result in the following files:
```ls -lh comparison/isec```

>``total 257K``<br>
>``-rw-r--r-- 1 bo1vsx bo 31K Feb 18 04:05 0000.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo 268 Feb 18 04:05 0000.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo 41K Feb 18 04:05 0001.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo 254 Feb 18 04:05 0001.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo 76K Feb 18 04:05 0002.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo 277 Feb 18 04:05 0002.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo 85K Feb 18 04:05 0003.bcf``<br>
>``-rw-r--r-- 1 bo1vsx bo 272 Feb 18 04:05 0003.bcf.csi``<br>
>``-rw-r--r-- 1 bo1vsx bo 559 Feb 18 04:05 README.txt``<br>

The content of the BCF files is explained in the README.txt file:
* ``0000.bcf``: private to bcftools
* ``0001.bcf``: private to gatk
* ``0002.bcf``: bcftools SNPs shared with gatk
* ``0003.bcf``: gatk SNPs shared with bcftools

### Filtering
The steps above got us a a few hundreds of SNPs, but usually we want to apply some filters depending on the sort of downstream analyses that we intend to do. Filtering usually consist on establishing some sort of threshold for the SNP quality in the QUAL field or some of the variables encoded in the INFO field (e.g. depth, mapping quality, etc.), but more sophisticated filtering is possible based on genotype-specific (i.e. per-sample) variables (e.g. discarding genotypes below a certain quality or with too few reads). Although in principle it is possible to do soft-filtering of SNPs and simply tag them as `PASSED` or `FAILED` in the VCF/BCF file, in practice it is more useful to do hard-filtering and simply remove all those SNPs (and genotypes) that did not meet our thresholds. We are going to focus on the SNPs only, so let's extract only the snps from all the variants shared by bcftools and gatk:

```bash
mkdir filtering
bcftools view -v snps comparison/isec/0002.bcf -O b > filtering/snps.bcf
bcftools index filtering/snps.bcf
bcftools view -H filtering/snps.bcf | wc -l
```

#### SNP-based filtering
bcftools allows applying filters on many of its commands, but usually they are used with ``bcftools view`` or with ``bcftools filter``. Filtering can be done using information encoded in the QUAL or INFO fields, also allowing expression with multiple conditions and basic arithmetics (more details [here](https://samtools.github.io/bcftools/bcftools.html#expressions)). These are some examples:

##### Filter by SNP quality
An obvious filter is to exclude (-e) SNPs below a quality threshold:
```bash
bcftools view -e 'QUAL<20' -O b filtering/snps.bcf > filtering/snps.QS20.bcf
```
Filters can also be specified as includes (-i), the equivalent of the one above is:
```bash
bcftools view -i 'QUAL>=20' -O b filtering/snps.bcf > filtering/snps.QS20i.bcf
```
You can see they have the same number of SNPs:
```bash
bcftools view -H filtering/snps.QS20.bcf | wc -l
bcftools view -H filtering/snps.QS20i.bcf | wc -l
```
Comparing the contents show they are identical:
```bash
diff -s <(bcftools view -H filtering/snps.QS20.bcf) <(bcftools view -H filtering/snps.QS20i.bcf)
```
##### Filter by SNP depth 
An example to exclude SNPs with depth <100:
```bash
bcftools view -e 'INFO/DP<100' -O b filtering/snps.bcf > filtering/snps.DP100.bcf
# You can see this is quite an stringent filter:
bcftools view -H filtering/snps.DP100.bcf | wc -l
```
An important note here is to be aware that if the file is further processed so that only part of the individuals are used, the fields in INFO may not be updated and it would be then unreliable to filter out using any information from that field.

##### Filter by mapping quality
Another typical filter is to remove SNPs with low mapping quality (RMS MQ: root mean square of the mapping quality of reads across all samples). For examples, all SNPs with RMS MQ < 20 can be discarded with the following command:
```bash
bcftools view -e 'MQ<20' -O b filtering/snps.bcf > filtering/snps.MQ20.bcf
```
If you check the number of SNPs, you can see this filter doesn't remove any:
```bash
bcftools view -H filtering/snps.bcf | wc -l
bcftools view -H filtering/snps.MQ20.bcf | wc -l
```
This is because all SNPs have a MQ=60:
```bash
bcftools query -f '%MQ\n' filtering/snps.bcf | less -S
bcftools query -f '%MQ\n' filtering/snps.bcf | uniq
```
##### Filter by allele frequency
For population genetic analyses it is frequently advised to remove variants below a certain allele frequency, as these ones are difficult to tell apart from sequencing errors. For example, to exclude all SNPs with a minor allele frequency (MAF) below 5% we would run:
```bash
bcftools view -e 'MAF<0.05' -O b filtering/snps.bcf > filtering/snps.MAF005.bcf
```
##### Filter by sample coverage
Another typical situation is to want to exclude all SNPs for which only a small fracion of all the individuals have sequence data. That requires using the total number of alleles in the called genotypes (AN). In the case of a diploid species, this can be used to filter SNPs by number of individuals genotyped. For example, we can remove all SNPs genotyped at less than 13 individuals:
```bash
bcftools view -e 'AN/2<13' -O b filtering/snps.bcf > filtering/snps.SAMP13.bcf
```

#### Genotype-based filtering
Filtering using the genotype fields can allow for some more precise filtering. For example, in cases of high depth heterogeneity among samples, it may be better to filter out by median genotype depth than by total depth across all samples. An example removing all SNPs where the mean genotype depth is below 5:
```bash
bcftools view -e 'AVG(FMT/DP)<5' -O b filtering/snps.bcf > filtering/snps.MEANGTDP5.bcf
```
Sometimes it is reasonable to ignore genotype calls based on few reads. The following command remove all genotype calls (i.e. genotypes are substituted by `./.`) informed by less than 3 reads:
```bash
bcftools filter -S . -e 'FMT/DP<3' -O b filtering/snps.bcf > filtering/snps.NOGTDP3.bcf
```
To see the effect of this, we can calculate the number and fraction of individuals genotyped for each SNP:
```bash
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n" filtering/snps.bcf | \
awk '{SUM=0; N=0; for(i=5; i<=NF; i++){if ($i!="./.") SUM+=1}; AVG=SUM/NF; print $1,$2,$3,$4,SUM,AVG}' | less -S

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n" filtering/snps.NOGTDP3.bcf | \
awk '{SUM=0; N=0; for(i=5; i<=NF; i++){if ($i!="./.") SUM+=1}; AVG=SUM/NF; print $1,$2,$3,$4,SUM,AVG}' | less -S
```

#### Combining multiple filters 
Multiple filters can be combined in a single command using or piping several ones. For example, we can combine a few of the filters we have used above:
```bash
bcftools filter -S . -e 'FMT/DP<3' filtering/snps.bcf | \
bcftools view -e 'AVG(FMT/DP)<5 || MAF<0.05 || QUAL<30 || AN/2<13' -O b > filtering/snps.NOGTDP3.MEANGTDP5.MAF005.Q30.SAMP13.bcf
# This results in quite a dramatic reduction in the number of SNPs:
bcftools view -H filtering/snps.NOGTDP3.MEANGTDP5.MAF005.Q30.SAMP13.bcf | wc -l
```
__Final note:__ It is important to be careful with the order of the filters, as different combinations can result in different end results, especially if we subsample the individuals at some point. Also, some positions may end up being invariant or only variant with respect to the reference (i.e. private), or not having any individuals genotyped (you could have spotted some examples above...).

[Back to TOC](index.md)
