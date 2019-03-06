rm(list=ls())
the.genotypes <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/Cowpea.GAPIT.GenoFormat.Lipka.txt", header=T)
the.genotypes[1:6,1:6]
the.genotypes$Snp <- paste("V", the.genotypes$Snp, sep="")


geno <- the.genotypes[,-c(2:5)]
rownames(geno) <- as.vector(as.matrix(geno[,1]))
geno <- geno[,-1]
geno[1:6,1:6]
geno <- as.matrix(geno)
geno[which(is.na(geno))] <- 1

G <- t(geno-1)
library(rrBLUP)

# Read in phenotypic data
phdata <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.phenotypes.csv", header = T)
head(phdata)
colnames(phdata)[1] <- "Taxa"
dim(phdata)


geno.taxa <- data.frame(colnames(geno))
head(geno.taxa)
colnames(geno.taxa)[1] <- "Taxa"
com.tax1 <- merge(phdata, geno.taxa, by="Taxa")
head(com.tax1)

# match genotypes by common taxa
the.genotypes <- the.genotypes[, c(1:5, match(com.tax1$Taxa, colnames(the.genotypes)))]
G[1:6,1:6]
G2 <- G[match(com.tax1$Taxa, rownames(G)),]
G2[1:6,1:6]
# Impute missing data using rrBLUP function
impute=A.mat(G2,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed

pheno1 <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/Cowpea.BLUPs.csv", header=T)
pheno2 <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.phenotypes.csv", header=T)
colnames(pheno2)[1] <- "Taxa"
pheno <- merge(pheno1, pheno2, by="Taxa")
head(pheno)
names(pheno)
#pheno <- pheno[,c(1,5,6,13,14,17,18)]
head(pheno)


geno.taxa <- data.frame(rownames(Markers_impute))
colnames(geno.taxa)[1] <- "Taxa"
head(geno.taxa)
com.tax <- merge(geno.taxa, pheno)
head(com.tax)

Markers_impute.t <- t(Markers_impute)


geno.epi <- Markers_impute.t[,match(com.tax$Taxa, colnames(Markers_impute.t))]

geno.epi[1:6,1:6]
class(geno.epi)
dim(geno.epi)

map <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.map.txt", header=T)
head(map)
map$LG.Pos <- paste(map$LG, map$Position, sep="_")
map2 <- map[!duplicated(map$LG.Pos),]
head(map2)
map2$Marker <- paste("V", map2$Marker, sep="")
head(map2)
dim(map2)
geno.epi.copy <- geno.epi
geno.epi[1:6,1:6]
geno.epi <- geno.epi+1
geno.epi[1:6,1:6]
geno.epi[geno.epi==1] <- 0.5
geno.epi[geno.epi==2] <- 1
geno.epi[1:6,1:6]

geno.epi.unique <- geno.epi#[match(map2$Marker, rownames(geno.epi)),]

dim(geno.epi.unique)
geno.epi.unique[1:6,1:6]

setwd("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/")

for(i in 2:ncol(com.tax)){
  trt.epi <- com.tax[,c(1,i)]
  trt.epi <- trt.epi[!is.na(trt.epi[,2]),]
  
  geno.trt.epi <- geno.epi.unique[,match(trt.epi$Taxa, colnames(geno.epi.unique))]
  geno.trt.epi <- t(geno.trt.epi)
  geno.trt.epi <- data.frame(geno.trt.epi)
  geno.trt.epi$Marker <- rownames(geno.trt.epi)
  geno.trt.epi <- geno.trt.epi[,c(ncol(geno.trt.epi), 2:(ncol(geno.trt.epi)-1))]
  colnames(geno.trt.epi)[1] <- "<Marker>"
  #geno.trt.epi <- as.matrix(geno.trt.epi)
  
  geno.trt.epi.nu <- geno.epi[,match(trt.epi$Taxa, colnames(geno.epi))]
  geno.trt.epi.nu <- t(geno.trt.epi.nu)
  geno.trt.epi.nu <- data.frame(geno.trt.epi.nu)
  geno.trt.epi.nu$Marker <- rownames(geno.trt.epi.nu)
  geno.trt.epi.nu <- geno.trt.epi.nu[,c(ncol(geno.trt.epi.nu), 2:(ncol(geno.trt.epi.nu)-1))]
  colnames(geno.trt.epi.nu)[1] <- "<Marker>"
  write.table(geno.trt.epi, paste("geno_TASSEL.JL_", colnames(com.tax)[i], ".txt", sep=""), sep="\t", quote = F, row.names=F, col.names=T)
  write.table(trt.epi, paste("TASSEL_", colnames(com.tax)[i], ".txt", sep=""), sep="\t", quote = F, row.names=F, col.names=T)
  
  write.table(geno.trt.epi.nu, paste("geno_TASSEL_", colnames(com.tax)[i], ".txt", sep=""), sep="\t", quote = F, row.names=F, col.names=T)
  #write.table(trt.epi[,2], paste("TASSEL_", colnames(com.tax)[i], ".txt", sep=""), sep="\t", quote = F, row.names=F, col.names=T)
}
write.table(map2, "cowpea.TASSEL.map.uniq.txt", sep="\t", quote=F, row.names=F, col.names=T)


########
####################### Epistasis Markers PVE and Additive Effects
# For Epistasis Markers
setwd("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL")

FT_BLUP <- read.csv("FT_BLUP.TASSEL.EpiMapping.Result.csv", header=T)
head(FT_BLUP)
FT_BLUP.1 <- FT_BLUP[,c(3,5,1)]
head(FT_BLUP.1)
colnames(FT_BLUP.1) <- c("SNP", "pval", "Trait")
FT_BLUP.2 <- FT_BLUP[,c(4,5,1)]
head(FT_BLUP.2)
colnames(FT_BLUP.2) <- c("SNP", "pval", "Trait")

FTFILD <- read.csv("FTFILD.TASSEL.EpiMapping.Result.csv", header=T)
head(FTFILD)
FTFILD.1 <- FTFILD[,c(3,5,1)]
colnames(FTFILD.1) <- c("SNP", "pval", "Trait")
FTFILD.2 <- FTFILD[,c(4,5,1)]
colnames(FTFILD.2) <- c("SNP", "pval", "Trait")
head(FTFILD.2)

FTRILD <- read.csv("FTRILD.TASSEL.EpiMapping.Result.csv", header=T)
head(FTRILD)
FTRILD.1 <- FTRILD[,c(3,5,1)]
colnames(FTRILD.1) <- c("SNP", "pval", "Trait")
FTRILD.2 <- FTRILD[,c(4,5,1)]
colnames(FTRILD.2) <- c("SNP", "pval", "Trait")
head(FTRILD.2)


FTFISD <- read.csv("FTFISD.TASSEL.EpiMapping.Result.csv", header=T)
head(FTFISD)
FTFISD.1 <- FTFISD[,c(3,5,1)]
colnames(FTFISD.1) <- c("SNP", "pval", "Trait")
FTFISD.2 <- FTFISD[,c(4,5,1)]
colnames(FTFISD.2) <- c("SNP", "pval", "Trait")
head(FTFISD.2)

FTRISD <- read.csv("FTRISD.TASSEL.EpiMapping.Result.csv", header=T)
head(FTRISD)
FTRISD.1 <- FTRISD[,c(3,5,1)]
colnames(FTRISD.1) <- c("SNP", "pval", "Trait")
FTRISD.2 <- FTRISD[,c(4,5,1)]
colnames(FTRISD.2) <- c("SNP", "pval", "Trait")
head(FTRISD.2)


MFISD <- read.csv("MFISD.TASSEL.EpiMapping.Result.csv", header=T)
head(MFISD)
MFISD.1 <- MFISD[,c(3,5,1)]
colnames(MFISD.1) <- c("SNP", "pval", "Trait")
MFISD.2 <- MFISD[,c(4,5,1)]
colnames(MFISD.2) <- c("SNP", "pval", "Trait")
head(MFISD.2)


MRISD <- read.csv("MRISD.TASSEL.EpiMapping.Result.csv", header=T)
head(MRISD)
MRISD.1 <- MRISD[,c(3,5,1)]
colnames(MRISD.1) <- c("SNP", "pval", "Trait")
MRISD.2 <- MRISD[,c(4,5,1)]
colnames(MRISD.2) <- c("SNP", "pval", "Trait")
head(MRISD.2)


MT_BLUP <- read.csv("MT_BLUP.TASSEL.EpiMapping.Result.csv", header=T)
head(MT_BLUP)
MT_BLUP.1 <- MT_BLUP[,c(3,5,1)]
colnames(MT_BLUP.1) <- c("SNP", "pval", "Trait")
MT_BLUP.2 <- MT_BLUP[,c(4,5,1)]
colnames(MT_BLUP.2) <- c("SNP", "pval", "Trait")
head(MT_BLUP.2)

SSFISD <- read.csv("SSFISD.TASSEL.EpiMapping.Result.csv", header=T)
head(SSFISD)
SSFISD.1 <- SSFISD[,c(3,5,1)]
colnames(SSFISD.1) <- c("SNP", "pval", "Trait")
SSFISD.2 <- SSFISD[,c(4,5,1)]
colnames(SSFISD.2) <- c("SNP", "pval", "Trait")
head(SSFISD.2)


SSRISD <- read.csv("SSRISD.TASSEL.EpiMapping.Result.csv", header=T)
head(SSRISD)
SSRISD.1 <- SSRISD[,c(3,5,1)]
colnames(SSRISD.1) <- c("SNP", "pval", "Trait")
SSRISD.2 <- SSRISD[,c(4,5,1)]
colnames(SSRISD.2) <- c("SNP", "pval", "Trait")
head(SSRISD.2)

SSRISD <- read.csv("SSRISD.TASSEL.EpiMapping.Result.csv", header=T)
head(SSRISD)
SSRISD.1 <- SSRISD[,c(3,5,1)]
colnames(SSRISD.1) <- c("SNP", "pval", "Trait")
SSRISD.2 <- SSRISD[,c(4,5,1)]
colnames(SSRISD.2) <- c("SNP", "pval", "Trait")
head(SSRISD.2)

SS_BLUP <- read.csv("SS_BLUP.TASSEL.EpiMapping.Result.csv", header=T)
head(SS_BLUP)
SS_BLUP.1 <- SS_BLUP[,c(3,5,1)]
colnames(SS_BLUP.1) <- c("SNP", "pval", "Trait")
SS_BLUP.2 <- SS_BLUP[,c(4,5,1)]
colnames(SS_BLUP.2) <- c("SNP", "pval", "Trait")
head(SS_BLUP.2)


JL_RES <- rbind(FT_BLUP.1, FT_BLUP.2, FTFILD.1, FTFILD.2, FTRILD.1, FTRILD.2, FTFISD.1, FTFISD.2, FTRISD.1, FTRISD.2, MFISD.1, MFISD.2, MRISD.1, MRISD.2, MT_BLUP.1, MT_BLUP.2, SSFISD.1, SSFISD.2, SSRISD.1, SSRISD.2, SS_BLUP.1, SS_BLUP.2)
JL_RES$Trait <- as.character(JL_RES$Trait)
JL_RES$SNP <- as.character(JL_RES$SNP)
head(JL_RES)


the.genotypes <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/Cowpea.GAPIT.GenoFormat.Lipka.txt", header=T)
the.genotypes[1:6,1:6]
the.genotypes$Snp <- paste("V", the.genotypes$Snp, sep="")


geno <- the.genotypes[,-c(2:5)]
rownames(geno) <- as.vector(as.matrix(geno[,1]))
geno <- geno[,-1]
geno[1:6,1:6]
geno <- as.matrix(geno)
geno[which(is.na(geno))] <- 1

G <- t(geno-1)
library(rrBLUP)

# Read in phenotypic data
phdata <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.phenotypes.csv", header = T)
head(phdata)
colnames(phdata)[1] <- "Taxa"
dim(phdata)


geno.taxa <- data.frame(colnames(geno))
head(geno.taxa)
colnames(geno.taxa)[1] <- "Taxa"
com.tax1 <- merge(phdata, geno.taxa, by="Taxa")
head(com.tax1)

# match genotypes by common taxa
the.genotypes <- the.genotypes[, c(1:5, match(com.tax1$Taxa, colnames(the.genotypes)))]
G[1:6,1:6]
G2 <- G[match(com.tax1$Taxa, rownames(G)),]
G2[1:6,1:6]
# Impute missing data using rrBLUP function
impute=A.mat(G2,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed

pheno1 <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/Cowpea.BLUPs.csv", header=T)
pheno2 <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.phenotypes.csv", header=T)
colnames(pheno2)[1] <- "Taxa"
pheno <- merge(pheno1, pheno2, by="Taxa")
head(pheno)
names(pheno)
pheno <- pheno[,-c(3,4,11,12,15,16)]
head(pheno)


geno.taxa <- data.frame(rownames(Markers_impute))
colnames(geno.taxa)[1] <- "Taxa"
head(geno.taxa)
com.tax <- merge(geno.taxa, pheno)
head(com.tax)
#colnames(com.tax)[2] <- "FTFILD"
Phenotypes.FT <- com.tax
# Impute missing data using rrBLUP function
head(com.tax)

phenames <- names(com.tax[,-1])

Markers_impute2 <- Markers_impute+1
Markers_impute2[1:6,1:6]
Markers_impute2 <- Markers_impute2[match(com.tax$Taxa, rownames(Markers_impute2)),]

setwd("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/")

for (l in 1:length(phenames))
{
  print(paste("-------------- Trait being analysed: ", phenames[l], "!!!!!!!!!!!---------------", sep = ""))
  
  ExplVar200Best <- JL_RES[which(JL_RES$Trait==phenames[l]),]
  
  bSNP<-Markers_impute2[,as.character(ExplVar200Best$SNP)]
  
  phdata <- data.frame(Phenotypes.FT[,1], Phenotypes.FT[,phenames[l]])
  colnames(phdata)[2] <- phenames[l]
  colnames(phdata)[1] <- "Taxa"
  #sP<-as.data.frame(phdata[,phenames[l]])
  sP<-phdata
  rownames(sP) <- sP$Taxa
  da<-as.data.frame(cbind(sP, bSNP))
  
  trait_QTL_Pheno <- da
  write.table(t(data.frame(c("Trait","QTL", "Additive Effect", "PVE"))), paste("Cowpea.Epistasis.QTL.Effects_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
  #APV is the Among population variance in accordance to Wurschum et al. 2011 Heredity
  for(i in 3:ncol(trait_QTL_Pheno)){
    
    snp <- colnames(trait_QTL_Pheno)[i]  
    print(paste("-------------- Trait being analysed: ", phenames[l], "SNP: ", snp, "!!!!!!!!!!!---------------", sep = ""))
    
    trait_QTL_Pheno_2 <- trait_QTL_Pheno[,c(1,2,i)]
    
    AA_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==2),]
    AA <- mean(AA_class[,2], na.rm=T)
    BB_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==0),]
    BB <- mean(BB_class[,2], na.rm=T)
    QTL_effect <- (AA-BB)/2
    
    #formula.single <- as.formula(paste("Cd_comb ~ ",paste(as.character(topSNP$SNP), collapse=" + "), sep=" "))
    trait_QTL_Pheno_2$QTL <- trait_QTL_Pheno_2[,3]
    #QTL <- colnames(trait_QTL_Pheno_2[3])
    fin.anova <- lm(trait_QTL_Pheno_2[,phenames[l]] ~  QTL, data=trait_QTL_Pheno_2, na.action = na.omit)
    fin.sum <- summary(fin.anova)
    
    QVar <- round((fin.sum$adj.r.squared)*100, digits=2)#Phenotypes.FT[,phenames[l]]
    
    print(paste("-------------- PVE For SNP: ", snp, "; Trait: ", phenames[l], " == ", QVar, "%  !!!!!!!!!!!---------------", sep = ""))
    write.table(t(data.frame(c(phenames[l],colnames(trait_QTL_Pheno[i]), round(abs(QTL_effect[1]), 1), QVar[1]))), paste("Cowpea.Epistasis.QTL.Effects_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
  }
  
  
}

# Dr. Lipka's model did a great job of identifying epistatic QTL
FT.BLUP.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_FLT_BLUP_QTL.txt", header=T)
FTFILD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_FTFILD_QTL.txt", header=T)
FTRILD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_FTRILD_QTL.txt", header=T)
FTFISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_FTFISD_QTL.txt", header=T)
FTRISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_FTRISD_QTL.txt", header=T)
MT.BLUP.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_MAT_BLUP_QTL.txt", header=T)
MFISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_MFISD_QTL.txt", header=T)
MRISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_MRISD_QTL.txt", header=T)
SS.BLUP.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_SS_BLUP_QTL.txt", header=T)
SSFISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_SSFISD_QTL.txt", header=T)
SSRISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.Epistasis.QTL.Effects_SSRISD_QTL.txt", header=T)

EpiQTL <- rbind(FT.BLUP.PVE, FTFILD.PVE, FTRILD.PVE, FTFISD.PVE, FTRISD.PVE, MT.BLUP.PVE, MFISD.PVE, MRISD.PVE, SS.BLUP.PVE, SSFISD.PVE, SSRISD.PVE)
head(EpiQTL)

EpiQTL.F <- rbind(FT.BLUP.PVE, FTFILD.PVE, FTRILD.PVE, FTFISD.PVE, FTRISD.PVE)
unique(EpiQTL.F$QTL)
EpiQTL.M <- rbind(MT.BLUP.PVE, MFISD.PVE, MRISD.PVE)
unique(EpiQTL.M$QTL)
EpiQTL.S <- rbind(SS.BLUP.PVE, SSFISD.PVE, SSRISD.PVE)
unique(EpiQTL.S$QTL)

Markers_impute2[1:6,1:6]
Markers_impute3 <- t(Markers_impute2)
Markers_impute3[1:6,1:6]
source("http://evachan.org/calc_snp_stats.R")

QTL.summary <- calc_snp_stats(Markers_impute3)
head(QTL.summary)
QTL.summary$QTL <- rownames(QTL.summary)
QTL.summary <- QTL.summary[,c(14,5,6)]
#hist(QTL.summary$maf)
summary(QTL.summary)

QTL.Info.Epi <- merge(EpiQTL, QTL.summary, by="QTL")
QTL.Info.Epi$maf <- round(QTL.Info.Epi$maf, 2)
head(QTL.Info.Epi)
dim(QTL.Info.Epi)
dim(JL_RES)
QTL.Info.Epi <- QTL.Info.Epi[order(QTL.Info.Epi$Trait, QTL.Info.Epi$PVE, decreasing = T),]
head(QTL.Info.Epi)
write.table(QTL.Info.Epi, "TASSEL.Epi.Summary.PVE.MAF.csv", sep=",", quote=F, row.names = F, col.names = T)

QTL.Info.Epi <- QTL.Info.Epi[,-5]
head(QTL.Info.Epi)

tables.vec <- c("FT_BLUP.TASSEL.EpiMapping.Result.csv", "FTFILD.TASSEL.EpiMapping.Result.csv", "FTRILD.TASSEL.EpiMapping.Result.csv",
                "FTFISD.TASSEL.EpiMapping.Result.csv", "FTRISD.TASSEL.EpiMapping.Result.csv", "MFISD.TASSEL.EpiMapping.Result.csv",
                "MRISD.TASSEL.EpiMapping.Result.csv", "MT_BLUP.TASSEL.EpiMapping.Result.csv", "SSFISD.TASSEL.EpiMapping.Result.csv",
                "SSRISD.TASSEL.EpiMapping.Result.csv", "SSRISD.TASSEL.EpiMapping.Result.csv", "SS_BLUP.TASSEL.EpiMapping.Result.csv")

for(j in 1:length(tables.vec)){
      
  trt.epis <- read.csv(tables.vec[j], header=T)
  head(trt.epis)
  trt.epis.1 <- trt.epis[,c(3,5,1,7,8)]
  head(trt.epis.1)
  colnames(trt.epis.1)[1] <- "QTL"
  nam.trt <- unique(as.character(trt.epis.1$Trait))
  QTL.Info.Epi$Trait <- as.character(QTL.Info.Epi$Trait)
  QTL.Info.Epi.trt <- QTL.Info.Epi[which(QTL.Info.Epi$Trait==nam.trt),]
  trt.epis.1 <- merge(trt.epis.1, QTL.Info.Epi.trt, by="QTL")
  head(trt.epis.1)
  trt.epis.1.2 <- trt.epis.1[match(trt.epis$Marker1, trt.epis.1$QTL),]
  trt.epis.1.2 <- trt.epis.1.2[,-6]
  colnames(trt.epis.1.2) <- c("Marker1", "P.Val", "Trait", "LG.1", "Position.1", "Additive.Effect.1", "PVE.1", "MAF.1")
  
  trt.epis.2 <- trt.epis[,c(4,5,1,10,11)]
  head(trt.epis.2)
  colnames(trt.epis.2)[1] <- "QTL"
  nam.trt2 <- unique(as.character(trt.epis.2$Trait))
  QTL.Info.Epi$Trait <- as.character(QTL.Info.Epi$Trait)
  QTL.Info.Epi.trt2 <- QTL.Info.Epi[which(QTL.Info.Epi$Trait==nam.trt2),]
  trt.epis.2 <- merge(trt.epis.2, QTL.Info.Epi.trt2, by="QTL")
  head(trt.epis.2)
  trt.epis.2.2 <- trt.epis.2[match(trt.epis$Marker2, trt.epis.2$QTL),]
  trt.epis.2.2 <- trt.epis.2.2[,-c(2,3,6)]
  names(trt.epis.2.2) <- c("Marker2", "LG.2", "Position.2", "Additive.Effect.1", "PVE.1", "MAF.1")
  trt.epis.info <- cbind(trt.epis.1.2, trt.epis.2.2)
  write.table(trt.epis.info, paste("Final.TASSEL.Epistasis.Result", nam.trt[1], "csv", sep="."), sep=",", quote=F, row.names = F, col.names=T)
  
}

head(QTL.Info.Epi)

map <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.map.txt", header=T)
head(map)
map$LG.Pos <- paste(map$LG, map$Position, sep="_")
map2 <- map[!duplicated(map$LG.Pos),]
head(map2)
map2$Marker <- paste("V", map2$Marker, sep="")
head(map2)
rownames(map2) <- NULL
dim(QTL.Info.Epi)
map4 <- map2
names(map4)[1] <- "QTL"

QTL.Info.Epi.f <- QTL.Info.Epi[grep("FLT", QTL.Info.Epi$Trait),]
QTL.Info.Epi.f <- QTL.Info.Epi[grep("FLT", QTL.Info.Epi$Trait),]
QTL.Info.Epi.s <- QTL.Info.Epi[grep("SS", QTL.Info.Epi$Trait),]
QTL.Info.Epi.m <- QTL.Info.Epi[grep("M", QTL.Info.Epi$Trait),]

QTL.Info.Epi.map <- merge(map4, QTL.Info.Epi, by="QTL")
QTL.Info.Epi.map <- QTL.Info.Epi.map[order(QTL.Info.Epi.map$Trait, QTL.Info.Epi.map$PVE, decreasing = T),]
head(QTL.Info.Epi.map, 15)

Links.info <- QTL.Info.Epi.map[,-c(1,4)]
Links.info.trt <- unique(as.character(Links.info$Trait))
Links.info2 <- Links.info
Links.info2$end <- Links.info2$Position+1.5

colnames(Links.info2)[1:2] <- c("chr", "start")  

for(m in 1:length(Links.info.trt)){
    Links.info.t <- Links.info2[which(Links.info2$Trait==Links.info.trt[m]),]
    maf.t <- Links.info.t[,c(1,2,7,6)]
    maf.t$color <- rep("a", nrow(maf.t))
    colnames(maf.t)[4] <- "value"
    write.table(maf.t, paste("Final.TASSEL.Epistasis.MAF", Links.info.trt[m], "csv", sep="."), sep=",", quote=F, row.names = F, col.names=T)
    
    
    adt.t <- Links.info.t[,c(1,2,7,4)]
    adt.t$color <- rep("b", nrow(adt.t))
    colnames(adt.t)[4] <- "value"
    write.table(adt.t, paste("Final.TASSEL.Epistasis.ADT", Links.info.trt[m], "csv", sep="."), sep=",", quote=F, row.names = F, col.names=T)
    
    
    pve.t <- Links.info.t[,c(1,2,7,5)]
    pve.t$color <- rep("c", nrow(pve.t))
    colnames(pve.t)[4] <- "value"
    write.table(pve.t, paste("Final.TASSEL.Epistasis.PVE", Links.info.trt[m], "csv", sep="."), sep=",", quote=F, row.names = F, col.names=T)
    
    
  }

######## Epistasis QTL Summary
FT_BLUP.Final <- read.csv("Final.TASSEL.Epistasis.Result.FLT_BLUP.csv", header=T)
FTFILD.Final <- read.csv("Final.TASSEL.Epistasis.Result.FTFILD.csv", header=T)
FTRILD.Final <- read.csv("Final.TASSEL.Epistasis.Result.FTRILD.csv", header = T)
FTFISD.Final <- read.csv("Final.TASSEL.Epistasis.Result.FTFISD.csv", header=T)
FTRISD.Final <- read.csv("Final.TASSEL.Epistasis.Result.FTRISD.csv", header=T)
MT_BLUP.Final <- read.csv("Final.TASSEL.Epistasis.Result.MAT_BLUP.csv", header = T)
MFISD.Final <- read.csv("Final.TASSEL.Epistasis.Result.MFISD.csv", header=T)
MRISD.Final <- read.csv("Final.TASSEL.Epistasis.Result.MRISD.csv", header=T)
SS_BLUP.Final <- read.csv("Final.TASSEL.Epistasis.Result.SS_BLUP.csv", header = T)
SSFISD.Final <- read.csv("Final.TASSEL.Epistasis.Result.SSFISD.csv", header=T)
SSRISD.Final <- read.csv("Final.TASSEL.Epistasis.Result.SSRISD.csv", header=T)

E_QTL <- rbind(FT_BLUP.Final, FTFILD.Final, FTRILD.Final, FTFISD.Final, FTRISD.Final, MT_BLUP.Final, MFISD.Final, MRISD.Final,
               SS_BLUP.Final, SSFISD.Final, SSRISD.Final)

#E_QTL <- read.csv("SPAEML.Epistatic.QTL.Summary.csv", header=T)
head(E_QTL)
E_QTL$QTL1 <- paste(E_QTL$LG.1, round(E_QTL$Position.1, 2), sep=":")
E_QTL$QTL1 <- paste("qVu", E_QTL$QTL1, sep="")
head(E_QTL)

head(E_QTL)
E_QTL$QTL2 <- paste(E_QTL$LG.2, round(E_QTL$Position.2, 2), sep=":")
E_QTL$QTL2 <- paste("qVu", E_QTL$QTL2, sep="")
head(E_QTL)
names(E_QTL)
E_QTL <- E_QTL[,c(1,15,2:8,9,16,10:14)]
head(E_QTL)
write.table(E_QTL, "SPAEML.Epistatic.QTL.Pos.Summary.New.02.09.2019.csv", sep=",", quote = F, row.names = F, col.names = T)


############################################################################
/homes/omo/tassel-5-standalone/./run_pipeline.pl -Xms249g -Xmx250g -fork1 -importGuess 
/homes/omo/public_html/Collaborations/Cowpea/TASSEL/geno_TASSEL_FTFILD.txt -fork2 -t 
/homes/omo/public_html/Collaborations/Cowpea/TASSEL/TASSEL_FTFILD.txt -combine3 -input1 -input2 -intersect 
-StepwiseOLSModelFitterPlugin -enter 0.001 -exit 0.001 -nestMarkers false -nperm 100 -endPlugin -export 
/homes/omo/public_html/Collaborations/Cowpea/TASSEL/FTFILD_stepwisePerm100.txt -runfork1 -runfork2


