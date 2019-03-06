############## Summarize Epistasis QTLs
rm(list=ls())
setwd("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL")
map3 <- read.table("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/cowpea.SIS.map.uniq.txt", header=T)
head(map3)
map3$SN <- c(1:nrow(map3))

bonferroni_05.epi <- -log10(1.556178e-06) #0.05/32130 #5.807941
bonferroni_05.epi2 <- -log10(1.556178e-05)

filenames= c("FT_BLUP.TASSEL.EpiResult.csv", "FTFILD.TASSEL.EpiResults.csv", "FTFISD.TASSEL.EpiResults.csv",
             "FTRILD.TASSEL.EpiResults.csv", "FTRISD.TASSEL.EpiResults.csv", "MT_BLUP.TASSEL.EpiResults.csv",
             "MRISD.TASSEL.EpiResults.csv", "MFISD.TASSEL.EpiResults.csv", "SS_BLUP.TASSEL.EpiResults.csv",
             "SSFISD.TASSEL.EpiResults.csv", "SSRISD.TASSEL.EpiResuts.csv")

traitnames= c("FT_BLUP", "FTFILD", "FTFISD","FTRILD", "FTRISD", "MT_BLUP","MRISD", "MFISD", "SS_BLUP","SSFISD", "SSRISD")

for(f in 1:length(filenames)){
trt <- read.table(filenames[f], sep=",", header=T)
trt <- trt[-1,]
trt$Marker1 <- sub('\\+.*', '', trt$Name)
trt$Marker2 <- sub('.*\\+', '', trt$Name)

trt <- trt[order(trt$pr.F),]
trt$p.log10 <- -log10(trt$pr.F)
trt <- trt[which(trt$p.log10>=6),]
trt <- trt[,c(1,2,15,16,10)]

trt.EpiMarker1 <- as.vector(trt$Marker1)
trt.Epi.info1 <- map3[which(map3$Marker %in%trt.EpiMarker1),c(1:3)]
head(trt.Epi.info1)
rownames(trt.Epi.info1) <- NULL
colnames(trt.Epi.info1) <- c("Marker1", "LG.1", "Position.1")
head(trt.Epi.info1)
trt.Epi.info1 <- trt.Epi.info1[match(trt$Marker1, trt.Epi.info1$Marker),]
#trt.Epi.info1.uniq <- trt.Epi.info1[!duplicated(trt.Epi.info1$LG.Pos.1),]
dim(trt.Epi.info1)


trt.EpiMarker2 <- as.vector(trt$Marker2)
trt.Epi.info2 <- map3[which(map3$Marker %in%trt.EpiMarker2),c(1:3)]
trt.Epi.info2.2 <- trt.Epi.info2[match(trt$Marker2, trt.Epi.info2$Marker),]
rownames(trt.Epi.info2.2) <- NULL
colnames(trt.Epi.info2.2) <- c("Marker2", "LG.2", "Position.2")
trt.Epi.info2.2

trt.Epi.Result <- cbind(trt.Epi.info1, trt.Epi.info2.2)
rownames(trt) <- NULL
trt.Epi.Result2 <- cbind(trt, trt.Epi.Result)
trt.Epi.Result2
trt.Epi.Result3 <- trt.Epi.Result[,c(2,3,5,6)]
colnames(trt.Epi.Result3) <- c("chr1", "start1", "chr2", "start2")
trt.Epi.Result3$end1 <- trt.Epi.Result3$start1+2
trt.Epi.Result3$end2 <- trt.Epi.Result3$start2+2
trt.Epi.Result3 <- trt.Epi.Result3[,c(1,2,5,3,4,6)]
write.table(trt.Epi.Result2, paste(traitnames[f], "TASSEL.EpiMapping.Result","csv", sep="."), sep=",", quote=F, row.names = F, col.names = T)
write.table(trt.Epi.Result3, paste(traitnames[f], "TASSEL.Links.Circos","csv", sep="."), sep=",", quote=F, row.names = F, col.names = T)
trt=NULL
}

head(map3)
write.table(map3, "Cowpea.UniqueGeneticMap.Epistasis.csv", sep=",", quote=F, row.names = F, col.names = T)


## SSRISD
SSRISD.inter <- read.table("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/MAT.SS.Epistasis/SSRISD.PATOWAS.Inter.uniq.txt", sep=",", header=F)
head(SSRISD.inter)
names(SSRISD.inter) <- c("Marker1", "Marker2", "Pval")

SSRISD.inter$logPval <- -log10(SSRISD.inter$Pval)
head(SSRISD.inter)

SSRISD.sig.bonf <- SSRISD.inter[which(SSRISD.inter$logPval >= bonferroni_05.epi),]
write.table(SSRISD.sig.bonf, "SSRISD.sig.bonf.csv", sep=",", quote=F, row.names = F, col.names = T)

# This was done manually
SSRISD.EpiMarker1 <- as.vector(SSRISD.sig.bonf$Marker1)
SSRISD.Epi.info1 <- map3[which(map3$SN %in%SSRISD.EpiMarker1),]
rownames(SSRISD.Epi.info1) <- NULL
colnames(SSRISD.Epi.info1) <- c("Marker1", "LG.1", "Position.1", "LG.Pos.1", "SN.1")
head(SSRISD.Epi.info1)

SSRISD.EpiMarker2 <- as.vector(SSRISD.sig.bonf$Marker2)
SSRISD.Epi.info2 <- map3[which(map3$SN %in%SSRISD.EpiMarker2),]
SSRISD.Epi.info2.1 <- SSRISD.Epi.info2[-c(1,2),]
SSRISD.Epi.info2.2 <- SSRISD.Epi.info2[match(SSRISD.sig.bonf$Marker2, SSRISD.Epi.info2$SN),]
rownames(SSRISD.Epi.info2.2) <- NULL
colnames(SSRISD.Epi.info2.2) <- c("Marker2", "LG.2", "Position.2", "LG.Pos.2", "SN.2")

SSRISD.Epi.Result <- cbind(SSRISD.Epi.info1, SSRISD.Epi.info2.2)
rownames(SSRISD.sig.bonf) <- NULL
SSRISD.Epi.Result2 <- cbind(SSRISD.sig.bonf, SSRISD.Epi.Result)
SSRISD.Epi.Result2
write.table(SSRISD.Epi.Result2, "SSRISD.EpiMapping.Result.csv", sep=",", quote=F, row.names = F, col.names = T)

# SS_BLUP
SS_BLUP.inter <- read.table("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/MAT.SS.Epistasis/MAT_BLUP.PATOWAS.Inter.uniq.txt", sep=",", header=F)
head(SS_BLUP.inter)
names(SS_BLUP.inter) <- c("Marker1", "Marker2", "Pval")

SS_BLUP.inter$logPval <- -log10(SS_BLUP.inter$Pval)
head(SS_BLUP.inter)

SS_BLUP.sig.bonf <- SS_BLUP.inter[which(SS_BLUP.inter$logPval >= bonferroni_05.epi),]
SS_BLUP.EpiMarker1 <- as.vector(SS_BLUP.sig.bonf$Marker1)
SS_BLUP.Epi.info1 <- map3[which(map3$SN %in%SS_BLUP.EpiMarker1),]
rownames(SS_BLUP.Epi.info1) <- NULL
colnames(SS_BLUP.Epi.info1) <- c("Marker1", "LG.1", "Position.1", "LG.Pos.1", "SN.1")

SS_BLUP.EpiMarker2 <- as.vector(SS_BLUP.sig.bonf$Marker2)
SS_BLUP.Epi.info2 <- map3[which(map3$SN %in%SS_BLUP.EpiMarker2),]
SS_BLUP.Epi.info2.2 <- SS_BLUP.Epi.info2[match(SS_BLUP.sig.bonf$Marker2, SS_BLUP.Epi.info2$SN),]
rownames(SS_BLUP.Epi.info2.2) <- NULL
colnames(SS_BLUP.Epi.info2.2) <- c("Marker2", "LG.2", "Position.2", "LG.Pos.2", "SN.2")

SS_BLUP.Epi.Result <- cbind(SS_BLUP.Epi.info1, SS_BLUP.Epi.info2.2)
rownames(SS_BLUP.sig.bonf) <- NULL
SS_BLUP.Epi.Result2 <- cbind(SS_BLUP.sig.bonf, SS_BLUP.Epi.Result)
SS_BLUP.Epi.Result2
write.table(SS_BLUP.Epi.Result2, "SS_BLUP.EpiMapping.Result.csv", sep=",", quote=F, row.names = F, col.names = T)

####################### Epistasis Markers PVE and Additive Effects
# For Epistasis Markers
trt <- read.csv("trt.EpiMapping.Result.csv", header=T)
head(trt)
trt$Trait <- rep("trt", nrow(trt))
trt.1 <- trt[,c(5,3,15)]
colnames(trt.1) <- c("SNP", "pval", "Trait")
head(trt.1)
trt.2 <- trt[,c(10,3,15)]
colnames(trt.2) <- c("SNP", "pval", "Trait")
head(trt.2)

SSRISD <- read.csv("SSRISD.EpiMapping.Result.csv", header=T)
head(SSRISD)
SSRISD$Trait <- rep("SSRISD", nrow(SSRISD))
SSRISD.1 <- SSRISD[,c(5,3,15)]
head(SSRISD.1)
colnames(SSRISD.1) <- c("SNP", "pval", "Trait")
SSRISD.2 <- SSRISD[,c(10,3,15)]
head(SSRISD.2)
colnames(SSRISD.2) <- c("SNP", "pval", "Trait")

SS_BLUP <- read.csv("SS_BLUP.EpiMapping.Result.csv", header=T)
head(SS_BLUP)
SS_BLUP$Trait <- rep("SS_BLUP", nrow(SS_BLUP))
SS_BLUP.1 <- SS_BLUP[,c(5,3,15)]
colnames(SS_BLUP.1) <- c("SNP", "pval", "Trait")
SS_BLUP.2 <- SS_BLUP[,c(10,3,15)]
colnames(SS_BLUP.2) <- c("SNP", "pval", "Trait")
head(SS_BLUP.2)

JL_RES <- rbind(trt.1, trt.2, SSRISD.1, SSRISD.2, SS_BLUP.1, SS_BLUP.2)
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
pheno <- pheno[,c(1,6,17,18)]
head(pheno)


geno.taxa <- data.frame(rownames(Markers_impute))
colnames(geno.taxa)[1] <- "Taxa"
head(geno.taxa)
com.tax <- merge(geno.taxa, pheno)
head(com.tax)
colnames(com.tax)[2] <- "SS_BLUP"
Phenotypes.FT <- com.tax
# Impute missing data using rrBLUP function
head(com.tax)

phenames <- names(com.tax[,-1])

Markers_impute2 <- Markers_impute+1
Markers_impute2[1:6,1:6]
Markers_impute2 <- Markers_impute2[match(com.tax$Taxa, rownames(Markers_impute2)),]

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
  write.table(t(data.frame(c("QTL", "Additive Effect", "PVE"))), paste("Cowpea.Epistasis.QTL.Effects_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
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
    write.table(t(data.frame(c(colnames(trait_QTL_Pheno[i]), round(abs(QTL_effect[1]), 1), QVar[1]))), paste("Cowpea.Epistasis.QTL.Effects_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
  }
  
  
}

trt.PVE.epi <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/MAT.SS.Epistasis/Cowpea.Epistasis.QTL.Effects_trt_QTL.txt", header=T)
trt.PVE.epi$Trait <- rep("trt", nrow(trt.PVE.epi))

SSRISD.PVE.epi <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/MAT.SS.Epistasis/Cowpea.Epistasis.QTL.Effects_SSRISD_QTL.txt", header=T)
SSRISD.PVE.epi$Trait <- rep("SSRISD", nrow(SSRISD.PVE.epi))

SS_BLUP.PVE.epi <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/MAT.SS.Epistasis/Cowpea.Epistasis.QTL.Effects_SS_BLUP_QTL.txt", header=T)
SS_BLUP.PVE.epi$Trait <- rep("SS_BLUP", nrow(SS_BLUP.PVE.epi))
SS.QTL.Epi <- rbind(trt.PVE.epi, SSRISD.PVE.epi, SS_BLUP.PVE.epi)
head(SS.QTL.Epi)

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

SS.Info.Epi <- merge(SS.QTL.Epi, QTL.summary, by="QTL")
SS.Info.Epi$maf <- round(SS.Info.Epi$maf, 2)
write.table(SS.Info.Epi, "SS.Epi.Summary.PVE.MAF.csv", sep=",", quote=F, row.names = F, col.names = T)








