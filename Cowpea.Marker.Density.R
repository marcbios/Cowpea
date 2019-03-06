rm(list=ls())

setwd("/homes/omo/public_html/Collaborations/Cowpea/Results/")
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) 
library("scatterplot3d")
library(kernlab)
library(BGLR)

source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/mlmm.r")
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/mlmm_cof.r")
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/plot_mlmm.r")
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sorghum/CoincidenceIndex.R")
#install.packages("https://github.com/Gregor-Mendel-Institute/mlmm/files/1356516/emma_1.1.2.tar.gz", repos = NULL)
library(emma) 

#Load your hapmap data and use GAPT to convert it to numeric format
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
phdata <- read.csv("http://people.beocat.ksu.edu/~omo/Collaborations/Cowpea/Cowpea.BLUPs.csv", header = T)
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


nsnps = 3 # number of SNPs to select from GWAS             #  enter heritability of the trait here
s <- 0.2  #  enter selection index here (proportion to select)                 
Iter=2500
burn=1000

#traits=1 # number of traits
cycles=1 # number of cross validation runs 100

Rx.Traits <- NULL
Ro.Traits <- NULL
Rk.Traits <- NULL
Sv.Trts <- NULL

dens.vec <- c(0.20,0.40,0.60,0.80)
seed.vec1 <- c(101:(100+cycles))
seed.vec2 <- c(201:(200+cycles))



for (l in 1:length(phenames)){
  print(paste("-------------- Trait being analysed: ", phenames[l], " !!!!!!!!!!!---------------", sep = ""))
 
   for(d in 1:length(dens.vec)){
     print(paste("-------------- Trait being analysed: ", phenames[l], "; Density: ",dens.vec[d]," !!!!!!!!!!!---------------", sep = ""))
        Rx.accuracy=matrix(nrow=cycles,ncol=4)
        Ro.accuracy=matrix(nrow=cycles,ncol=4)
        Rk.accuracy=matrix(nrow=cycles,ncol=4)
        Sv.accuracy=matrix(nrow=cycles,ncol=4)
       
  
  
        Pheno=as.matrix(com.tax[,phenames[l]])
        rownames(Pheno) <- com.tax$Taxa
        
        GEN <- Markers_impute
        GENO <- Markers_impute + 1
        M <-tcrossprod(GENO)/ncol(GENO)
        X <- GENO
        
          for(r in 1:cycles)
          {
            print(paste("-------------- Trait being analysed: ", phenames[l], "; Density: ",dens.vec[d]," --cycle-- ", r, " !!!!!!!!!!!---------------", sep = ""))
            #define training and validation set
            #set.seed(seed.vec1[r])
            train=as.matrix(sample(1:nrow(Pheno),dens.vec[d]*nrow(Pheno)), replace=F) # number of genotypes selected out of total, here 80% --> chnage 0.8 if wanted
            valid<-setdiff(1:nrow(Pheno),train)
            
            #set.seed(seed.vec2[r])
            valid2 <- sample(valid, 193, replace = F)
            # taining set
            Pheno_train=Pheno[train,] # phenotypes
            m_train1=G[train,]
            
            
            #validation set 
            Pheno_valid=Pheno[valid2,] # phenotype
            m_valid1=G[valid2,]
           
            Y = Pheno_train
            Y=Y[!is.na(Y)]
            genot.trait=m_train1[match(names(Y),rownames(m_train1)),]
            
            k2 <-Kmat[match(names(Y), rownames(Kmat)),]
            k22 <- as.matrix(k2)
            k3 <-k22[,match(names(Y), colnames(k22))]
            
            
            print(paste("-------------- Running MLMM for: ", phenames[l], " --cycle-- ", r, " !!!!!!!!!!!---------------", sep = ""))
            
            mygwas_trait<-mlmm(Y=Y,X=genot.trait, K=k3,2,10)
            res.Trait=mygwas_trait$opt_mbonf$out
            res.Trait2 <- res.Trait[order(res.Trait$pval, decreasing = F),]
            head(res.Trait2)
            
            Gwas.output<- res.Trait2
            
            these.markers <- as.character(Gwas.output[1:nsnps,1])
            
            best <- these.markers
            # taining set
            m_train=m_train1[, !as.character(colnames(m_train1)) %in% best] # get marker data without fixed effect markers
            m_ft=m_train1[, as.character(colnames(m_train1)) %in% best] # markers as fixed
            m_fix_t <- cbind(as.matrix(rep(1, length(Pheno_train))), m_ft) # make matrix with mean (1) and fixed markers
            colnames(m_fix_t)=c("mean",best)
            
            # validation set 
            m_valid=m_valid1[, !colnames(m_valid1) %in% best] # markers without fixed markers
            m_fv=m_valid1[, colnames(m_valid1) %in% best] # markers as fixed
            m_fix_v <- cbind(as.matrix(rep(1, length(Pheno_valid))), m_fv) # make matrix with mean (1) and fixed markers
            colnames(m_fix_v)=c("mean",best)
            
            
            #FxRRBLUP
            print(paste("--------- Trait being analysed: ", phenames[l], " --cycle-- ", r, " ; GS Model: FxRRBLUP"," !!!!!!!----------", sep = ""))
            library(rrBLUP)
            Yt <-Pheno_train
            
            rrMod.Rx<-mixed.solve(Yt, X=m_fix_t, Z=m_train, K=NULL, SE=F, return.Hinv=F)
            mEff.Rx<-rrMod.Rx$u
            e.Rx= as.matrix(mEff.Rx)
            
            predYv.Rx = m_valid%*%e.Rx
            
            predYr.Rx = predYv.Rx[,1]+ (m_fix_v %*% rrMod.Rx$beta)
            
            Y_valid=Pheno_valid
            Rx.accuracy[r,1] <- cor(predYr.Rx,Y_valid,use="complete")
            the.fitted.model.Rx <- lm(Y_valid ~ predYr.Rx)
            the.coefficients.Rx <- c(the.fitted.model.Rx$coefficients[1], the.fitted.model.Rx$coefficients[2])
            Rx.accuracy[r,2] <- the.coefficients.Rx[1]
            Rx.accuracy[r,3] <- the.coefficients.Rx[2]
            
            these.observed.and.predicted.phenotypic.values.Rx <- data.frame(rownames(predYr.Rx), Y_valid, predYr.Rx)
            rownames(these.observed.and.predicted.phenotypic.values.Rx) <- NULL
            colnames(these.observed.and.predicted.phenotypic.values.Rx) <- c("Taxa", "Observed.Value", "Predicted.Value")
            x.p.Rx=these.observed.and.predicted.phenotypic.values.Rx[,c(1,3)]
            y.o.Rx=these.observed.and.predicted.phenotypic.values.Rx[,c(1,2)]
            Rx.accuracy[r,4] <- round(CI(x.p.Rx,y.o.Rx,s=s,top=T),2)
            
            
            
            #Ordinary RRBLUP
            print(paste("--------- Trait being analysed: ", phenames[l], "; Density: ",dens.vec[d], " --cycle-- ", r, " ; GS Model: RRBLUP"," !!!!!!!----------", sep = ""))
            rrMod.Ro<-mixed.solve(Yt, X=NULL, Z=m_train1, K=NULL, SE=F, return.Hinv=F)
            mEff.Ro<-rrMod.Ro$u
            e.Ro= as.matrix(mEff.Ro)
            
            predYv.Ro = m_valid1%*%e.Ro
            
            predYr.Ro = predYv.Ro[,1]+ rrMod.Ro$beta
            
            Y_valid=Pheno_valid
            Ro.accuracy[r,1] <- cor(predYr.Ro,Y_valid,use="complete")
            the.fitted.model.Ro <- lm(Y_valid ~ predYr.Ro)
            the.coefficients.Ro <- c(the.fitted.model.Ro$coefficients[1], the.fitted.model.Ro$coefficients[2])
            Ro.accuracy[r,2] <- the.coefficients.Ro[1]
            Ro.accuracy[r,3] <- the.coefficients.Ro[2]
            
            these.observed.and.predicted.phenotypic.values.Ro <- data.frame(names(predYr.Ro), Y_valid, predYr.Ro)
            rownames(these.observed.and.predicted.phenotypic.values.Ro) <- NULL
            colnames(these.observed.and.predicted.phenotypic.values.Ro) <- c("Taxa", "Observed.Value", "Predicted.Value")
            x.p.Ro=these.observed.and.predicted.phenotypic.values.Ro[,c(1,3)]
            y.o.Ro=these.observed.and.predicted.phenotypic.values.Ro[,c(1,2)]
            Ro.accuracy[r,4] <- round(CI(x.p.Ro,y.o.Ro,s=s,top=T),2)
            
            
            
            #RKHS 
            print(paste("--------- Trait being analysed: ", phenames[l], "; Density: ",dens.vec[d], " --cycle-- ", r, " ; GS Model: RKHS"," !!!!!!!----------", sep = ""))
            ETA<-list(list(K=M,model='RKHS')) 
            fm.RK<-BGLR(y=y.trn,ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
            Rk.accuracy[r,1] <- cor(fm.RK$yHat[valid2], ww, use="complete")
            the.fitted.model.RK <- lm(ww ~ fm.RK$yHat[valid2])
            the.coefficients.RK <- c(the.fitted.model.RK$coefficients[1], the.fitted.model.RK$coefficients[2])
            Rk.accuracy[r,2] <- the.coefficients.RK[1]
            Rk.accuracy[r,3] <- the.coefficients.RK[2]
            
            these.observed.and.predicted.phenotypic.values.Rk <- data.frame(names(Y_valid), Y_valid, fm.RK$yHat[valid2])
            rownames(these.observed.and.predicted.phenotypic.values.Rk) <- NULL
            colnames(these.observed.and.predicted.phenotypic.values.Rk) <- c("Taxa", "Observed.Value", "Predicted.Value")
            x.p.Rk=these.observed.and.predicted.phenotypic.values.Rk[,c(1,3)]
            y.o.Rk=these.observed.and.predicted.phenotypic.values.Rk[,c(1,2)]
            Rk.accuracy[r,4] <- round(CI(x.p.Rk,y.o.Rk,s=s,top=T),2)
            
            
            #SVR
            print(paste("--------- Trait being analysed: ", phenames[l], "; Density: ",dens.vec[d], " --cycle-- ", r, " ; GS Model: SVR"," !!!!!!!----------", sep = ""))
            svp_w <- ksvm(m_train1, Yt, type="eps-svr", kernel = "rbfdot")
            yhat <- predict(svp_w, m_valid1)
            rownames(yhat) <- names(Y_valid)
            Sv.accuracy[r,1] <- cor(yhat, Y_valid,use="complete")
            the.fitted.model.Sv <- lm(Y_valid ~ yhat)
            the.coefficients.Sv <- c(the.fitted.model.Sv$coefficients[1], the.fitted.model.Sv$coefficients[2])
            Sv.accuracy[r,2] <- the.coefficients.Sv[1]
            Sv.accuracy[r,3] <- the.coefficients.Sv[2]
            
            
            these.observed.and.predicted.phenotypic.values.Sv <- data.frame(names(Y_valid), Y_valid, yhat)
            rownames(these.observed.and.predicted.phenotypic.values.Sv) <- NULL
            colnames(these.observed.and.predicted.phenotypic.values.Sv) <- c("Taxa", "Observed.Value", "Predicted.Value")
            x.p.Sv=these.observed.and.predicted.phenotypic.values.Sv[,c(1,3)]
            y.o.Sv=these.observed.and.predicted.phenotypic.values.Sv[,c(1,2)]
            Sv.accuracy[r,4] <- round(CI(x.p.Sv,y.o.Sv,s=s,top=T),2)
            
            
          } # End Cycles
          
          print(paste("--------- Summarizing Result for trait: ", phenames[l] , "; Density: ",dens.vec[d]," !!!!!!!----------", sep = ""))
          Rx.accuracy <- data.frame(Rx.accuracy)
          Rx.accuracy$Trt <- rep(phenames[l], nrow(Rx.accuracy))
          Rx.Traits <- rbind(Rx.Traits, Rx.accuracy)
          
          Ro.accuracy <- data.frame(Ro.accuracy)
          Ro.accuracy$Density <- rep(dens.vec[d], nrow(Ro.accuracy))
          Ro.Density <- rbind(Ro.Density, Ro.accuracy)
          
          Rk.accuracy <- data.frame(Rk.accuracy)
          Rk.accuracy$Density <- rep(dens.vec[d], nrow(Rk.accuracy))
          Rk.Density <- rbind(Rk.Density, Rk.accuracy)
          
          Sv.accuracy <- data.frame(Sv.accuracy)
          Sv.accuracy$Density <- rep(dens.vec[d], nrow(Sv.accuracy))
          Sv.Density <- rbind(Sv.Density, Sv.accuracy)
          
          
        } # Density
  print(paste("--------- Summarizing Result for trait: ", phenames[l]," !!!!!!!----------", sep = ""))
  
        Rx.Density <- data.frame(Rx.Density)
        Rx.Density$Traits <- rep(phenames[l], nrow(Rx.Density))
        Rx.Traits <- rbind(Rx.Traits, Rx.Density)
  
        Ro.Density <- data.frame(Ro.Density)
        Ro.Density$Traits <- rep(phenames[l], nrow(Ro.Density))
        Ro.Traits <- rbind(Ro.Traits, Ro.Density)
        
        
        Rk.Density <- data.frame(Rk.Density)
        Rk.Density$Traits <- rep(phenames[l], nrow(Rk.Density))
        Rk.Traits <- rbind(Rk.Traits, Rk.Density)
        
        Sv.Density <- data.frame(Sv.Density)
        Sv.Density$Traits <- rep(phenames[l], nrow(Sv.Density))
        Sv.Traits <- rbind(Sv.Traits, Sv.Density)

        Rx.Density = NULL
        Ro.Density = NULL
        Rk.Density = NULL
        Sv.Density = NULL
        

} # End of Traits


save.image("GS.Results.Density.rda")
#Rx.Traits$Model <- rep("wRRBLUP", nrow(Rx.Traits))
Ro.Traits$Model <- rep("RRBLUP", nrow(Ro.Traits))
Rk.Traits$Model <- rep("RKHS", nrow(Rk.Traits))
Sv.Traits$Model <- rep("SVR", nrow(Sv.Traits))

GS.Results <- rbind(Ro.Traits, Rk.Traits, Sv.Traits) #Rx.Traits, 
colnames(GS.Results)[c(1:4)] <- c("P.Accuracy", "Intercept", "Slope", "C.Index")
GS.Results <- GS.Results[order(GS.Results$Traits, GS.Results$Model),]
write.table(GS.Results, "GS.Results.New.csv", sep=",", quote=F, row.names=F, col.names=T)

GS.Results.PH.LP <- GS.Results[which(GS.Results$Trt=="PH.LP"),]
GS.Results.PH.HP <- GS.Results[which(GS.Results$Trt=="PH.HP"),]
GS.Results.DTFL.LP <- GS.Results[which(GS.Results$Trt=="DTFL.LP"),]
GS.Results.DTFL.HP <- GS.Results[which(GS.Results$Trt=="DTFL.HP"),]
GS.Results.GYLD.LP <- GS.Results[which(GS.Results$Trt=="GYLD.LP"),]
GS.Results.GYLD.HP <- GS.Results[which(GS.Results$Trt=="GYLD.HP"),]
GS.Results.PANL.LP <- GS.Results[which(GS.Results$Trt=="PANL.LP"),]
GS.Results.PANL.HP <- GS.Results[which(GS.Results$Trt=="PANL.HP"),]

GS.Pheno.Names <- c("PH.LP", "PH.HP", "DTFL.LP", "DTFL.HP", "GYLD.LP", "GYLD.HP", "PANL.LP", "PANL.HP")
pdf("GScomparison.Accuracy.new.pdf", 15, 45)
par(mfrow=c(4,2))
for(k in 1:length(GS.Pheno.Names)){
  GS.Results.Trt <- GS.Results[which(GS.Results$Trt==GS.Pheno.Names[k]),]
  boxplot(GS.Results.Trt$P.Accuracy ~  GS.Results.Trt$Model,col=c("light blue", "coral1", "forestgreen", "green", "limegreen", "gold1", "darkviolet", "khaki2"),
          ylab="Prediction accuracy", main=GS.Pheno.Names[k], cex.main=2, cex.lab=2, cex.axis=1.5)
}
dev.off()


pdf("GScomparison.CoincidenceIndex.pdf", 10, 20)
par(mfrow=c(3,2))
for(k in 1:length(GS.Pheno.Names)){
  GS.Results.Trt <- GS.Results[which(GS.Results$Trt==GS.Pheno.Names[k]),]
  boxplot(GS.Results.Trt$C.Index ~  GS.Results.Trt$Model,col=c("light blue", "coral1", "green", "gold1", "darkviolet", "khaki2"),
          ylab="Coincidence Index", main=GS.Pheno.Names[k], cex.main=2, cex.lab=2, cex.axis=1.5)
}
dev.off()


pdf("GScomparison.Intercept.Bias.pdf", 10, 20)
par(mfrow=c(3,2))
for(k in 1:length(GS.Pheno.Names)){
  GS.Results.Trt <- GS.Results[which(GS.Results$Trt==GS.Pheno.Names[k]),]
  boxplot(GS.Results.Trt$Intercept ~  GS.Results.Trt$Model,col=c("light blue", "coral1", "green", "gold1", "darkviolet", "khaki2"),
          ylab="Intercept", main=GS.Pheno.Names[k], cex.main=2, cex.lab=2, cex.axis=1.5)
}
dev.off()


pdf("GScomparison.Slope.Bias.pdf", 10, 20)
par(mfrow=c(3,2))
for(k in 1:length(GS.Pheno.Names)){
  GS.Results.Trt <- GS.Results[which(GS.Results$Trt==GS.Pheno.Names[k]),]
  boxplot(GS.Results.Trt$Slope ~  GS.Results.Trt$Model,col=c("light blue", "coral1", "green", "gold1", "darkviolet", "khaki2"),
          ylab="Slope", main=GS.Pheno.Names[k], cex.main=2, cex.lab=2, cex.axis=1.5)
}
dev.off()

