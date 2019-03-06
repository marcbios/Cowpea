source("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Functions/geno2allel.R")
df <- read.csv("cowpea.genome.mxn.nomap.csv", header=T)
df[1:6,1:6]
rownames(df) <- df$SNP.Name
df <- df[,-1]
df[1:6,1:6]
df <- as.matrix(df)
geno.numeric <- geno_to_allelecnt(df, ref=NULL)
dim(geno.numeric)
geno.numeric[1:6,1:6]
geno.numerict <- t(geno.numeric)
geno.numerict[1:6,1:6]
geno.numerict <- data.frame(geno.numerict)
geno.numerict[1:6,1:6]
write.table(geno.numerict, "cowpea.geno.txt", sep="\t", quote=F, row.names=T, col.names = T)

maf.cw <- calc_snp_stats(geno.numeric)
head(maf.cw)
summary(maf.cw$maf)# minimum minor allele frequency is 0.05254

geno.numeric[1:6,1:6]
class(geno.numeric)
geno.epi <- as.matrix(geno.numeric)
geno.epi[1:10,1:6]
geno.epi[which(geno.epi=="NA")] <- 1
geno.epi[which(is.na(geno.epi))] <- 1
class(geno.epi) <- "numeric"
geno.epi[1:10,1:6]
#geno.epi2 <- geno.epi + 1
geno.epi[geno.epi==0] <- -1
geno.epi[geno.epi==1] <- 0 
geno.epi[geno.epi==2] <- 1

write.table(geno.epi, "cowpea.PATOWAS.geno.csv", sep=",", quote=F, row.names = F, col.names = F)


#pheno <- read.csv("Cowpea.BLUPs.csv", header=T)
pheno <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/cowpea.phenotypes.csv", header=T)
head(pheno)
pheno <- pheno[,c(1,2,3)]
head(pheno)
colnames(pheno)[1] <- "Taxa"

geno.taxa <- data.frame(rownames(geno.numerict))
colnames(geno.taxa)[1] <- "Taxa"
head(geno.taxa)
com.tax <- merge(geno.taxa, pheno)
head(com.tax)

geno.numerict2 <- geno.numerict[match(rownames(geno.numerict), com.tax$Taxa),]
geno.epi <- geno.epi[,match(com.tax$Taxa, colnames(geno.epi))]
geno.epi2 <- geno.numeric[,match(com.tax$Taxa, colnames(geno.numeric))]

geno.epi[1:6,1:6]
class(geno.epi)

geno.epi[1:6,1:6]
geno.epi2[1:6,1:6]

geno.epi2[which(geno.epi2==2)] <- 0
for(i in 2:ncol(com.tax)){
  trt.epi <- com.tax[,c(1,i)]
  trt.epi <- trt.epi[!is.na(trt.epi[,2]),]
  geno.trt.epi2 <- geno.epi2[,match(trt.epi$Taxa, colnames(geno.epi2))]
  geno.trt.epi <- geno.epi[,match(trt.epi$Taxa, colnames(geno.epi))]
  write.table(geno.trt.epi2, paste("geno_pepis.dom_", colnames(com.tax)[i], ".csv", sep=""), sep=",", quote = F, row.names=F, col.names=F)
  write.table(geno.trt.epi, paste("geno_pepis.adt_", colnames(com.tax)[i], ".csv", sep=""), sep=",", quote = F, row.names=F, col.names=F)
  write.table(trt.epi[,2], paste("pepis_", colnames(com.tax)[i], ".csv", sep=""), sep=",", quote = F, row.names=F, col.names=F)
}


the.genotypes <- read.delim("http://people.beocat.ksu.edu/~omo/Collaborations/Cowpea/Cowpea.GAPIT.GenoFormat.Lipka.txt", header=T)
the.genotypes[1:6,1:6]
the.genotypes$Snp <- paste("V", the.genotypes$Snp, sep="")
snp_info<-the.genotypes[c(1,3:4)] 
colnames(snp_info)<-c("SNP","Chr","Pos")
head(snp_info)
map <- read.delim("cowpea.map.txt", header=T)
head(map)
map$LG.Pos <- paste(map$LG, map$Position, sep="_")
map2 <- map[!duplicated(map$LG.Pos),]
head(map2)

geno.epi[1:6,1:6]
geno.epi2[1:6,1:6]

geno.epi.unique <- geno.epi[match(map2$Marker, rownames(geno.epi)),]
geno.epi2.unique <- geno.epi2[match(map2$Marker, rownames(geno.epi2)),]

dim(geno.epi.unique)
dim(geno.epi2.unique)

for(i in 2:ncol(com.tax)){
  trt.epi <- com.tax[,c(1,i)]
  trt.epi <- trt.epi[!is.na(trt.epi[,2]),]
  #geno.trt.epi2 <- geno.epi2.unique[,match(trt.epi$Taxa, colnames(geno.epi2.unique))]
  geno.trt.epi <- geno.epi.unique[,match(trt.epi$Taxa, colnames(geno.epi.unique))]
  #write.table(geno.trt.epi2, paste("geno_pepis.dom.uniq_", colnames(com.tax)[i], ".csv", sep=""), sep=",", quote = F, row.names=F, col.names=F)
  write.table(geno.trt.epi, paste("geno_pepis.adt.uniq_", colnames(com.tax)[i], ".csv", sep=""), sep=",", quote = F, row.names=F, col.names=F)
  write.table(trt.epi[,2], paste("pepis_", colnames(com.tax)[i], ".csv", sep=""), sep=",", quote = F, row.names=F, col.names=F)
}


files.name <- c("YLD.BLUP.PEPIS_Result.txt", "FLT.BLUP.PEPIS_Result.txt", "HGT.BLUP.PEPIS_Result.txt", "TILLER.BLUP.PEPIS_Result.txt")

mat <- matrix(NA,length(files.name), 7)
mat <- data.frame(mat)
colnames(mat) <- c("Aditive", "Domiance", "Additive.Additive", "Additive.Domiance", "Domiance.Additive", "Domiance.Domiance", "Residual")
for(i in 1:length(files.name)){
  
  trait.epis.vc <- read.delim(files.name[i], header=T)
  for(j in 1:ncol(trait.epis.vc)){
    vc <- trait.epis.vc[1,j]
    vcp <- (vc/rowSums(trait.epis.vc))*100
    mat[i,j] <- round(vcp,1)
  }
  
}
mat.copy <- mat
mat$Trait <- c("YLD", "FLT", "HGT", "TILLER")
write.table(mat, "PEPIS.Rice.VC.txt", sep="\t", quote=F, row.names=F, col.names=T)

setwd("~/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results")


#files.name <- c("FLT.PATOWAS.VC.uniq.txt","GH.PATOWAS.VC.uniq.txt", "MAT.PATOWAS.VC.uniq.txt", "SS.PATOWAS.VC.uniq.txt", "YLD.PATOWAS.VC.uniq.txt")
files.name <- c("FTFILD.PATOWAS.VC.uniq.txt","FTRILD.PATOWAS.VC.uniq.txt", "FTFISD.PATOWAS.VC.uniq.txt", "FTRISD.PATOWAS.VC.uniq.txt", "FLT.PATOWAS.VC.uniq.txt")
#mat$Trait <- c("FTFILD", "FTRILD", "FTFISD", "FTRISD", "FT_BLUP")
mat <- matrix(NA,length(files.name), 3)
mat <- data.frame(mat)
for(i in 1:length(files.name)){
  
  trait.epis.vc <- read.table(files.name[i], sep=",", header=T)
  for(j in 1:ncol(trait.epis.vc)){
    vc <- trait.epis.vc[1,j]
    vcp <- (vc/rowSums(trait.epis.vc))*100
    mat[i,j] <- round(vcp,1)
  }
  
}
mat
names(mat) <- c("Additive", "Additive-Additive", "Residual")
mat
#mat.copy
#mat$Trait <- c("FTFILD", "FTFISD", "FTRILD", "FTRISD", "GHFISD", "GHRISD", "GYFISD", "GYRISD", "MFISD")
#mat$Trait <- c("FLT_BLUP","GH_BLUP", "MAT_BLUP", "SS_BLUP", "YLD_BLUP")
mat$Trait <- c("FTFILD", "FTRILD", "FTFISD", "FTRISD", "FT_BLUP")


write.table(mat, "PATOWAS.Cowpea.GenoUniq.FT.VC..txt", sep="\t", quote=F, row.names=F, col.names=T)
mat

############### Genome Position

get_genome_pos <- function( chr_num, position, buffer=10){
  # chr_num: vector with chromosome numbers
  # position: vector with positions, corresponding to each chr_num
  # buffer: space between chromosomes in bp
  # Chromosome lengths for sorghum:
  chr_length <- c(
    chromosome_1  = 86.3426,
    chromosome_2  = 72.5654,
    chromosome_3  = 132.6856,  
    chromosome_4  = 78.7016,
    chromosome_5  = 106.926,
    chromosome_6  = 80.6053,
    chromosome_7  = 104.9475,
    chromosome_8  = 78.1136,
    chromosome_9  = 97.1424,
    chromosome_10 = 70.8373,
    chromosome_11 = 70.6164
  )
  chr_length <- chr_length + buffer
  
  position + sapply( chr_num, function(x) sum(c(0,chr_length[-length(chr_length)])[ 1:x ]))
}

###############
FLT.MP <- read.table("FLT.PATOWAS.Main.uniq.txt", sep=",", header=F)
head(FLT.MP)
names(FLT.MP) <- c("SN", "P.value")
map2$SN <- c(1:nrow(map2))
FLT.MP2 <- merge(map2, FLT.MP, by="SN")
head(FLT.MP2)

Trait.MP <- FLT.MP2

Trait.MP$chr <- as.numeric(Trait.MP$LG )
Trait.MP$pos <- as.numeric(Trait.MP$Position)

Trait.MP <- Trait.MP[order(Trait.MP$chr, Trait.MP$pos),]
head(Trait.MP)

GH.MP <- read.table("GH.PATOWAS.Main.uniq.txt", sep=",", header=F)
head(GH.MP)
names(GH.MP) <- c("SN", "P.value")
map2$SN <- c(1:nrow(map2))
GH.MP2 <- merge(map2, GH.MP, by="SN")
head(GH.MP2)

GH.MP <- GH.MP2

GH.MP$chr <- as.numeric(GH.MP$LG )
GH.MP$pos <- as.numeric(GH.MP$Position)

GH.MP <- GH.MP[order(GH.MP$chr, GH.MP$pos),]
head(GH.MP)

GH.MP$genome_pos <- get_genome_pos(GH.MP$chr, GH.MP$pos)
head(GH.MP)


Trait.MP$genome_pos <- get_genome_pos(Trait.MP$chr, Trait.MP$pos)
head(Trait.MP)


bonferroni_05 <- -log10(0.05/nrow(Trait.MP))

mean_genome_pos <- tapply(Trait.MP$genome_pos, Trait.MP$chr, function(x) (min(x)+max(x))/2 )
vek<-as.numeric(Trait.MP$chr)%%2


mars <- c(1,7,2,0)
mars2 <- c(6,7,2,2)
man_pch=20

man_lwd=2
cexp2=2
cexs=2
col_vec_man <- c("black","gray45")

pdf("FLT.BLUP.PATOWAS.pdf", width = 12, height = 5)
par(mar=mars2)
  
plot( I(-log10(Trait.MP$P.value)) ~ Trait.MP$genome_pos, cex=cexp2, col=col_vec_man[as.factor(vek)], 
        ylab=expression(paste('-log'[10],(italic(p)))),  cex.lab=cexs, cex.main=cexs, cex.axis=cexs,
        xlab='Chromosome', xaxt='n', main='', pch=man_pch, frame.plot =F)

plot( I(-log10(GH.MP$P.value)) ~ GH.MP$genome_pos, cex=cexp2, col=col_vec_man[as.factor(vek)], 
      ylab=expression(paste('-log'[10],(italic(p)))),  cex.lab=cexs, cex.main=cexs, cex.axis=cexs,
      xlab='Chromosome', xaxt='n', main='', pch=man_pch, frame.plot =F)

#abline(h=bonferroni_05, lty='dashed', col="red", lwd=man_lwd)
axis(side=1, at=mean_genome_pos, labels=1:11, tick=F, cex.axis=cexs)

dev.off()


######################
map <- read.delim("cowpea.map.txt", header=T)
head(map)
map$LG.Pos <- paste(map$LG, map$Position, sep="_")
map2 <- map[!duplicated(map$LG.Pos),]
head(map2)
map2$SN <- c(1:nrow(map2))

bonferroni_05.epi <- -log10(1.556178e-06) #0.05/32130 #5.807941

#FTFILD.inter <- read.table("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/FTFILD.PATOWAS.Inter.uniq.txt", sep=",", header=F)
FTFILD.inter <- read.table("/homes/omo/public_html/Collaborations/Cowpea/Results/FTFILD.PATOWAS.Inter.2019.txt", sep=",", header=F)


head(FTFILD.inter)
names(FTFILD.inter) <- c("Marker1", "Marker2", "Pval")

FTFILD.inter$logPval <- -log10(FTFILD.inter$Pval)
head(FTFILD.inter)

FTFILD.sig.bonf <- FTFILD.inter[which(FTFILD.inter$logPval >= bonferroni_05.epi),]

FTFILD.sig.bonf <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results/FTFILD.sig.bonf.csv", header = T)
head(FTFILD.sig.bonf)
dim(FTFILD.sig.bonf)
FTFILD.EpiMarker1 <- as.vector(FTFILD.sig.bonf$Marker1)
FTFILD.Epi.info1 <- map2[which(map2$SN %in%FTFILD.EpiMarker1),]
dim(FTFILD.Epi.info1)
rownames(FTFILD.Epi.info1) <- NULL
colnames(FTFILD.Epi.info1) <- c("Marker1", "LG.1", "Position.1", "LG.Pos.1", "SN.1")
head(FTFILD.Epi.info1)
#FTFILD.Epi.info1.uniq <- FTFILD.Epi.info1[!duplicated(FTFILD.Epi.info1$LG.Pos.1),]
dim(FTFILD.Epi.info1)


FTFILD.EpiMarker2 <- as.vector(FTFILD.sig.bonf$Marker2)
FTFILD.Epi.info2 <- map2[which(map2$SN %in%FTFILD.EpiMarker2),]
FTFILD.Epi.info2.2 <- FTFILD.Epi.info2[match(FTFILD.sig.bonf$Marker2, FTFILD.Epi.info2$SN),]
rownames(FTFILD.Epi.info2.2) <- NULL
colnames(FTFILD.Epi.info2.2) <- c("Marker2", "LG.2", "Position.2", "LG.Pos.2", "SN.2")
FTFILD.Epi.info2.2

FTFILD.Epi.Result <- cbind(FTFILD.Epi.info1, FTFILD.Epi.info2.2)
rownames(FTFILD.sig.bonf) <- NULL
FTFILD.Epi.Result2 <- cbind(FTFILD.sig.bonf, FTFILD.Epi.Result)
FTFILD.Epi.Result2
write.table(FTFILD.Epi.Result2, "FTFILD.EpiMapping.Result.csv", sep=",", quote=F, row.names = F, col.names = T)
head(map2)
write.table(map2, "Cowpea.UniqueGeneticMap.Epistasis.csv", sep=",", quote=F, row.names = F, col.names = T)


## FTRILD
FTRILD.inter <- read.table("/homes/omo/public_html/Collaborations/Cowpea/Results/FTRILD.PATOWAS.Inter.2019.txt", sep=",", header=F)
head(FTRILD.inter)
names(FTRILD.inter) <- c("Marker1", "Marker2", "Pval")

FTRILD.inter$logPval <- -log10(FTRILD.inter$Pval)
head(FTRILD.inter)

FTRILD.sig.bonf <- FTRILD.inter[which(FTRILD.inter$logPval >= bonferroni_05.epi),]
FTRILD.EpiMarker1 <- as.vector(FTRILD.sig.bonf$Marker1)
FTRILD.Epi.info1 <- map2[which(map2$SN %in%FTRILD.EpiMarker1),]
rownames(FTRILD.Epi.info1) <- NULL
colnames(FTRILD.Epi.info1) <- c("Marker1", "LG.1", "Position.1", "LG.Pos.1", "SN.1")

FTRILD.EpiMarker2 <- as.vector(FTRILD.sig.bonf$Marker2)
FTRILD.Epi.info2 <- map2[which(map2$SN %in%FTRILD.EpiMarker2),]
FTRILD.Epi.info2.2 <- FTRILD.Epi.info2[match(FTRILD.sig.bonf$Marker2, FTRILD.Epi.info2$SN),]
rownames(FTRILD.Epi.info2.2) <- NULL
colnames(FTRILD.Epi.info2.2) <- c("Marker2", "LG.2", "Position.2", "LG.Pos.2", "SN.2")

FTRILD.Epi.Result <- cbind(FTRILD.Epi.info1, FTRILD.Epi.info2.2)
rownames(FTRILD.sig.bonf) <- NULL
FTRILD.Epi.Result2 <- cbind(FTRILD.sig.bonf, FTRILD.Epi.Result)
FTRILD.Epi.Result2
#write.table(FTRILD.Epi.Result2, "FTRILD.EpiMapping.Result.csv", sep=",", quote=F, row.names = F, col.names = T)

# FT_BLUP
FT_BLUP.inter <- read.table("/homes/omo/public_html/Collaborations/Cowpea/Results/FT_BLUP.PATOWAS.Inter.2019.txt", sep=",", header=F)
head(FT_BLUP.inter)
names(FT_BLUP.inter) <- c("Marker1", "Marker2", "Pval")

FT_BLUP.inter$logPval <- -log10(FT_BLUP.inter$Pval)
head(FT_BLUP.inter)

FT_BLUP.sig.bonf <- FT_BLUP.inter[which(FT_BLUP.inter$logPval >= bonferroni_05.epi),]
FT_BLUP.EpiMarker1 <- as.vector(FT_BLUP.sig.bonf$Marker1)
FT_BLUP.Epi.info1 <- map2[which(map2$SN %in%FT_BLUP.EpiMarker1),]
rownames(FT_BLUP.Epi.info1) <- NULL
colnames(FT_BLUP.Epi.info1) <- c("Marker1", "LG.1", "Position.1", "LG.Pos.1", "SN.1")

FT_BLUP.EpiMarker2 <- as.vector(FT_BLUP.sig.bonf$Marker2)
FT_BLUP.Epi.info2 <- map2[which(map2$SN %in%FT_BLUP.EpiMarker2),]
FT_BLUP.Epi.info2.2 <- FT_BLUP.Epi.info2[match(FT_BLUP.sig.bonf$Marker2, FT_BLUP.Epi.info2$SN),]
rownames(FT_BLUP.Epi.info2.2) <- NULL
colnames(FT_BLUP.Epi.info2.2) <- c("Marker2", "LG.2", "Position.2", "LG.Pos.2", "SN.2")

FT_BLUP.Epi.Result <- cbind(FT_BLUP.Epi.info1, FT_BLUP.Epi.info2.2)
rownames(FT_BLUP.sig.bonf) <- NULL
FT_BLUP.Epi.Result2 <- cbind(FT_BLUP.sig.bonf, FT_BLUP.Epi.Result)
FT_BLUP.Epi.Result2
write.table(FT_BLUP.Epi.Result2, "FT_BLUP.EpiMapping.Result.csv", sep=",", quote=F, row.names = F, col.names = T)







