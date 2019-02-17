## PCA of genoypes

We are going to carry out PCA using R. For that, we will need to convert our BCF file to mean genotype format (based on the BIMBAM format), where genotypes are encoded as mean genotypes. A mean genotype is a value between 0 and 2 that can be interpreted as the minor allele dosage: 0 is homozygous for the major allele, 1 is a heterozygote, and 2 is a homozygote for the minor allele. Thus, intermediate values reflect the uncertainty in genotypes. 

We are going to use a custom Perl script called ``bcf2bbgeno.pl`` to get such mean genotypes. This scripts (1) removes all the SNPs that have more than two alleles and (2) calculates empirical posterior genotype probabilities from the genotype likelihoods in the BCF file under the assumption that the population is in Hardy-Weinberg equilibrium (HWE). Specifically, the script uses inferred allele frequencies to set HWE priors:

p(AA) = p2; p(aa) = (1-p)2; p(Aa) = 2p(1-p)

being p the allele frequency of major/reference allele A. Genotype likelihoods are multiplied by these priors to obtain genotype posterior probabilities that are then encoded as mean genotypes and saved to a .bbgeno file.

You can get some info about how to run the Perl script:
```bash
mkdir scripts
# download script
wget https://raw.githubusercontent.com/visoca/popgenomworkshop-gwas_gemma/master/scripts/bcf2bbgeno.pl -O scripts/bcf2bbgeno.pl
# give execution permissions
chmod +x scripts/bcf2bbgeno.pl
# show help
scripts/bcf2bbgeno.pl -h
```

To calculate the genotype posterior probabilites and save them in mean genotype format, we would need to run a command like this one, including the flag to include a header with the individuals names (followed by compression to save some space, given that ```R``` can handle gzipped files without problems):
```bash
mkdir pca
scripts/bcf2bbgeno.pl -i filtering/snps.bcf -o pca/snps.bbgeno -p H-W -a
```
We will also need the information about the race and sex of the samples, which can be copied from a shared directory in ShARC:
```bash
cp
```
Then we will use R to do PCA with the mean genotypes and do some plotting. First, let's load the genotypes
```R
## load genotypes
genofile<-"pca/snps.bbgeno"
genotypes<-read.table(genofile,header=T, check.names=F)
```
Now we will calculate the genotype covariance matrix:
```R
## calculate N x N genotype covariance matrix
gmn<-apply(genotypes[,-(1:3)],1,mean)
nids<-ncol(genotypes)-3
nloci<-nrow(genotypes)
gmnmat<-matrix(gmn,nrow=nloci,ncol=nids)
gprime<-genotypes[,-(1:3)]-gmnmat

gcovarmat<-matrix(NA,nrow=nids,ncol=nids)
for(i in 1:nids){
 for(j in i:nids){
  if (i==j){
    gcovarmat[i,j]<-cov(gprime[,i],gprime[,j])
  }
  else{
    gcovarmat[i,j]<-cov(gprime[,i],gprime[,j])
    gcovarmat[j,i]<-gcovarmat[i,j]
  }
 }
}
```
And we run a PCA using the function `prcomp`:
```R
## pca on the genotype covariance matrix
pcg<-prcomp(x=gcovarmat,center=TRUE,scale=FALSE)
```
Extract the PCS:
```R
pcs<-pcg$x
rownames(pcs)<-colnames(genotypes[,-(1:3)])
pcs.col<-unlist(lapply (rownames(pcs), function(x)
as.character(col.table[col.table[,2]==populations[populations[,1]==as.character(x),][[2]],][[1]])))
```
```R

# plot each PC using different colours for race and different symbols for sex
race.col<-c("red","blue")
sex.symbol<-c(0,2)
#Plot PC1 vs PC2
plot(pcs[,1],pcs[, i+2], main="PC1 vs PC2", xlab="PC1", ylab="PC2"))
```

We can also plot the names of the samples to investigate potential outliers:
```R
# plot each PC with sample names
for (i in 1:(npcs-1)){
par(mfrow=c(1,1),oma=c(1,1,1,10), mar=c(5, 4, 4, 2) + 0.1)
plot(pcs[,i],pcs[,i+1],type="n", main="PCA using covariance matrix", xlab=paste("PC",i,sep=""), ylab=paste("PC",(i+1),sep=""))
text(pcs[,i],pcs[,i+1],labels=rownames(pcs),col=pcs.col,cex=0.5)
# plot legend to the right
par(mfrow=c(1,1), oma=c(0,0,0,1), mar=c(0,0,0,0), new=T)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
legend("right", legend=pops, lty=1, lwd=10, col=colours, bty="n", xpd=T)

```
