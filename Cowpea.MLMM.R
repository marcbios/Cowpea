
rm(list=ls())
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/mlmm.r")
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/mlmm_cof.r")
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/plot_mlmm.r")
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sorghum/CoincidenceIndex.R")
#install.packages("https://github.com/Gregor-Mendel-Institute/mlmm/files/1356516/emma_1.1.2.tar.gz", repos = NULL)
library(emma) 

# REad in Genetic data
the.genotypes <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/Cowpea.GAPIT.GenoFormat.Lipka.txt", header=T)
the.genotypes[1:6,1:6]
the.genotypes$Snp <- paste("V", the.genotypes$Snp, sep="")
snp_info<-the.genotypes[c(1,3:4)] 
colnames(snp_info)<-c("SNP","Chr","Pos")
head(snp_info)

geno <- the.genotypes[,-c(2:5)]
rownames(geno) <- as.vector(as.matrix(geno[,1]))
geno <- geno[,-1]
geno[1:6,1:6]
geno <- as.matrix(geno)
geno[which(is.na(geno))] <- 1
#geno[which(geno=="N")] <- 1
G <- t(geno-1)
library(rrBLUP)

# Read in phenotypic data
phdata <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.phenotypes.csv", header = T)
head(phdata)
colnames(phdata)[1] <- "Taxa"
dim(phdata)
phenames <- as.vector(colnames(phdata[,-1]))

geno.taxa <- data.frame(colnames(geno))
head(geno.taxa)
colnames(geno.taxa)[1] <- "Taxa"
com.tax <- merge(phdata, geno.taxa, by="Taxa")
head(com.tax)

# match genotypes by common taxa
the.genotypes <- the.genotypes[, c(1:5, match(com.tax$Taxa, colnames(the.genotypes)))]
G[1:6,1:6]
G2 <- G[match(com.tax$Taxa, rownames(G)),]

# Impute missing data using rrBLUP function
impute=A.mat(G2,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed
# Create Kinship matrix using rrBLUP
Kmat <- A.mat(Markers_impute, min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,shrink=FALSE,return.imputed=FALSE)

for(i in 2:ncol(com.tax)){
  
  pheno2 <- com.tax[,c(1,i)]
  trait.name <- names(pheno2[2])
  
  print(paste("-------------- Performing GWAS on trait ", trait.name, "!!!!!!!---------------", sep = "")) 

  Y = as.vector(pheno2[,2])
  names(Y) <- as.vector(as.matrix(pheno2$Taxa))
  Y=Y[!is.na(Y)]
  genot.trait=Markers_impute[match(names(Y),rownames(Markers_impute)),]
  
  k2 <-Kmat[match(names(Y), rownames(Kmat)),]
  k22 <- as.matrix(k2)
  k3 <-k22[,match(names(Y), colnames(k22))]
  
  #tryCatch({
    mygwas_trait<-mlmm(Y=Y,X=genot.trait, K=k3,2,10)
  #}, error=function(e){})
  #mygwas_trait<-mlmm(Y=Y,X=genot.trait, K=k3,2,10)
  res.Trait=mygwas_trait$opt_mbonf$out
  res.Trait2 <- res.Trait[order(res.Trait$pval, decreasing = F),]
  head(res.Trait2)
  write.table(res.Trait, paste("GWAS_mlmm_results", trait.name, "csv",sep="."), sep=",", quote=F, row.names=F, col.names=T)
  
  pdf(paste("opt_GWAS_mbonf_",trait.name,".pdf",sep=""),width=15,height=4,paper='special')
  plot_opt_GWAS(mygwas_trait,'mbonf',snp_info,0.1)
  dev.off()
  
  pdf(paste("qqplot_opt_GWAS_mbonf_",trait.name,".pdf",sep=""))
  qqplot_opt_GWAS(mygwas_trait,'mbonf')
  dev.off()

}






