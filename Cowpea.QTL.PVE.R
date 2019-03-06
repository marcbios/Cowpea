the.genotypes <- read.delim("http://people.beocat.ksu.edu/~omo/Collaborations/Cowpea/Cowpea.GAPIT.GenoFormat.Lipka.txt", header=T)
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
phdata <- read.csv("http://people.beocat.ksu.edu/~omo/Collaborations/Cowpea/cowpea.phenotypes.csv", header = T)
head(phdata)
colnames(phdata)[1] <- "Taxa"
dim(phdata)
phenames <- as.vector(colnames(phdata[,-1]))

geno.taxa <- data.frame(colnames(geno))
head(geno.taxa)
colnames(geno.taxa)[1] <- "Taxa"
com.tax <- merge(phdata, geno.taxa, by="Taxa")
head(com.tax)
com.tax.FT.SEnvs <- #com.tax[,c(1:5)]
com.tax.FT.BLUP <- #com.tax[,c(1,2)]
com.tax <- merge(com.tax.FT.SEnvs, com.tax.FT.BLUP, by="Taxa")
# match genotypes by common taxa
the.genotypes <- the.genotypes[, c(1:5, match(com.tax$Taxa, colnames(the.genotypes)))]
G[1:6,1:6]
G2 <- G[match(com.tax$Taxa, rownames(G)),]

Phenotypes.FT <- com.tax
# Impute missing data using rrBLUP function
impute=A.mat(G2,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed
Markers_impute[1:6,1:6]
Markers_impute <- Markers_impute
Markers_impute2 <- Markers_impute+1
Markers_impute2[1:6,1:6]

FTFILD <- read.csv("FTFILD.MLMM.csv", header=T)
FTRILD <- read.csv("FTRILD.MLMM.csv", header=T)
FTFISD <- read.csv("FTFISD.MLMM.csv", header=T)
FTRISD <- read.csv("FTRISD.MLMM.csv", header=T)
FT_BLUP <- read.csv("FT_BLUP.MLMM.csv", header=T)

JL_RES <- rbind(FTFILD, FTRILD, FTFISD, FTRISD, FT_BLUP)
JL_RES$Trait <- as.character(JL_RES$Trait)
JL_RES$SNP <- as.character(JL_RES$SNP)

phenames <- names(com.tax[,-1])
for (l in 1:length(phenames))
  #for(l in 1:3)
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
  write.table(t(data.frame(c("QTL", "Additive Effect", "PVE"))), paste("Cowpea.QTL.Effects_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
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
    write.table(t(data.frame(c(colnames(trait_QTL_Pheno[i]), round(abs(QTL_effect[1]), 1), QVar[1]))), paste("Cowpea.QTL.Effects_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
  }
  
  
}


FTFILD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/Cowpea.QTL.Effects_FTFILD_QTL.txt", header=T)
FTFILD.PVE$Trait <- rep("FTFILD", nrow(FTFILD.PVE))

FTRILD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/Cowpea.QTL.Effects_FTRILD_QTL.txt", header=T)
FTRILD.PVE$Trait <- rep("FTRILD", nrow(FTRILD.PVE))

FTFISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/Cowpea.QTL.Effects_FTFISD_QTL.txt", header=T)
FTFISD.PVE$Trait <- rep("FTFISD", nrow(FTFISD.PVE))

FTRISD.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/Cowpea.QTL.Effects_FTRISD_QTL.txt", header=T)
FTRISD.PVE$Trait <- rep("FTRISD", nrow(FTRISD.PVE))

FT_BLUP.PVE <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/Cowpea.QTL.Effects_FLT_BLUP_QTL.txt", header=T)
FT_BLUP.PVE$Trait <- rep("FT_BLUP", nrow(FT_BLUP.PVE))

FT.Eff <- rbind(FTFILD.PVE, FTRILD.PVE, FTFISD.PVE, FTRISD.PVE, FT_BLUP.PVE)


###################################################################################################################

# For Epistasis Markers
FTFILD <- read.csv("FTFILD.PATOWAS.csv", header=T)
FT_BLUP <- read.csv("FT_BLUP.PATOWAS.csv", header=T)

JL_RES <- rbind(FTFILD, FT_BLUP)
JL_RES$Trait <- as.character(JL_RES$Trait)
JL_RES$SNP <- as.character(JL_RES$SNP)

phenames <- names(com.tax[,c(2,6)])
for (l in 1:length(phenames))
  #for(l in 1:3)
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

FTFILD.PVE.epi <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/Cowpea.Epistasis.QTL.Effects_FTFILD_QTL.txt", header=T)
FTFILD.PVE.epi$Trait <- rep("FTFILD", nrow(FTFILD.PVE.epi))

FT_BLUP.PVE.epi <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/Cowpea.Epistasis.QTL.Effects_FLT_BLUP_QTL.txt", header=T)
FT_BLUP.PVE.epi$Trait <- rep("FT_BLUP", nrow(FT_BLUP.PVE.epi))
FT.QTL.Epi <- rbind(FTFILD.PVE.epi, FT_BLUP.PVE.epi)

############################### Minor Allele Frequency
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

FT.info <- merge(FT.Eff, QTL.summary, by="QTL")

head(FT.info)
FT.info$maf <- round(FT.info$maf, 2)
head(FT.info)

map <- read.delim("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.map.txt", header=T)
head(map)
map$Marker <- paste("V", map$Marker, sep="")
head(map)
map.q <- map
colnames(map.q)[1] <- "QTL"
FT.info2 <- merge(map.q, FT.info, by="QTL")
head(FT.info2)

write.table(FT.info2, "FT.PVE.ADE.MAF.MLMM.csv", sep=",", quote=F, row.names = F, col.names = T)

FT.Info.Epi <- merge(FT.QTL.Epi, QTL.summary, by="QTL")
FT.Info.Epi$maf <- round(FT.Info.Epi$maf, 2)
write.table(FT.Info.Epi, "FT.Epi.Summary.PVE.MAF.csv", sep=",", quote=F, row.names = F, col.names = T)




