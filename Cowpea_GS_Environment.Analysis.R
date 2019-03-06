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

FTFILD <- com.tax[,c(1,2)]
FTFILD$WStatus <- rep("FI", nrow(FTFILD))
FTFILD$PhotoPd <- rep("LD", nrow(FTFILD))
FTFILD$Env <- rep("FILD", nrow(FTFILD))
colnames(FTFILD)[2] <- "FLT"

FTRILD <- com.tax[,c(1,3)]
FTRILD$WStatus <- rep("RI", nrow(FTRILD))
FTRILD$PhotoPd <- rep("LD", nrow(FTRILD))
FTRILD$Env <- rep("RILD", nrow(FTRILD))
colnames(FTRILD)[2] <- "FLT"


FTFISD <- com.tax[,c(1,4)]
FTFISD$WStatus <- rep("FI", nrow(FTFISD))
FTFISD$PhotoPd <- rep("SD", nrow(FTFISD))
FTFISD$Env <- rep("FISD", nrow(FTFISD))
colnames(FTFISD)[2] <- "FLT"


FTRISD <- com.tax[,c(1,5)]
FTRISD$WStatus <- rep("RI", nrow(FTRISD))
FTRISD$PhotoPd <- rep("SD", nrow(FTRISD))
FTRISD$Env <- rep("RISD", nrow(FTRISD))
colnames(FTRISD)[2] <- "FLT"

FLT.Data <- rbind(FTFILD, FTRILD, FTFISD, FTRISD)
head(FLT.Data)
library(ggplot2)

# Flowering time reaction norm in cowpea
pdf("FLT.Rxn.Norm.2.pdf", 20,14)
ggplot(data=FLT.Data, aes(x=Env, y=FLT, group=Taxa)) +
  geom_line(color="wheat2", size=0.6)+
  geom_point(color="wheat4", size=2) +
  xlab("Environment") + ylab("Flowering Time") + labs(colour="Env")  + 
  theme_bw() + # remove grey background (because Tufte said so)
  theme(plot.title = element_text(size=size.title, face = "bold"), legend.title=element_text(size=size.leg,face="bold"),
        legend.text=element_text(colour="black", size=size.leg,face="bold"), 
        axis.title.x = element_text(face="bold", colour="black", size=size.axtit), 
        axis.title.y = element_text(face="bold", colour="black", size=size.axtit)
        , axis.text.x = element_text(color="black", size=size.axtxt, angle=angle.val),
        axis.text.y = element_text(color="black", size=size.axtxt), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

size.axtit=50
size.axtxt=40
angle.val=35
size.title=40
size.leg=15

pdf("FLT.Rxn.Norm.pdf", 20,14)
ggplot(data=FLT.Data, aes(x=Env, y=FLT, group=Taxa)) +
  geom_line(aes(color=FLT), size=0.8) +
  scale_color_gradient(low = "blue", high = "coral1") +
  xlab("Environment") + ylab("Flowering Time") + labs(colour="Env") +
  theme_bw() + # remove grey background (because Tufte said so)
  theme(plot.title = element_text(size=size.title, face = "bold"), legend.title=element_text(size=size.leg,face="bold"),
        legend.text=element_text(colour="black", size=size.leg,face="bold"), 
        axis.title.x = element_text(face="bold", colour="black", size=size.axtit), 
        axis.title.y = element_text(face="bold", colour="black", size=size.axtit)
        , axis.text.x = element_text(color="black", size=size.axtxt, angle=angle.val),
        axis.text.y = element_text(color="black", size=size.axtxt), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) #+ scale_fill_manual(values=group.colors)
dev.off()


pheno.b <- read.csv("http://people.beocat.ksu.edu/~omo/Collaborations/Cowpea/Cowpea.BLUPs.csv", header=T)
dim(pheno.b)
library(corrplot)
M<-cor(com.tax[,c(2:5)],use="complete")
head(round(M,2))
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(com.tax[,c(2:5)])
head(p.mat[, 1:5])
pdf("Cowpea.Env.FLT.Cor.pdf",8.2,11.6)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
dev.off()


M <- cor(pheno.b[,c(2:6)], use="complete")
p.mat <- cor.mtest(pheno.b[,c(2:6)])
head(p.mat[, 1:5])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
pdf("Cowpea.Pheno.Cor.pdf",8.2,11.6)
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
dev.off()

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
cycles=100 # number of cross validation runs 100

Rx.Trts <- NULL
Ro.Trts <- NULL
Rk.Trts <- NULL
Sv.Trts <- NULL

seed.vec1 <- c(101:(100+cycles))
Rx.accuracy=matrix(NA, length(phenames), length(phenames))
Ro.accuracy=matrix(NA, length(phenames), length(phenames))
Rk.accuracy=matrix(NA, length(phenames), length(phenames))
Sv.accuracy=matrix(NA, length(phenames), length(phenames))

Rx.slope=matrix(NA, length(phenames), length(phenames))
Ro.slope=matrix(NA, length(phenames), length(phenames))
Rk.slope=matrix(NA, length(phenames), length(phenames))
Sv.slope=matrix(NA, length(phenames), length(phenames))

Rx.intercept=matrix(NA, length(phenames), length(phenames))
Ro.intercept=matrix(NA, length(phenames), length(phenames))
Rk.intercept=matrix(NA, length(phenames), length(phenames))
Sv.intercept=matrix(NA, length(phenames), length(phenames))

Rx.coincidence=matrix(NA, length(phenames), length(phenames))
Ro.coincidence=matrix(NA, length(phenames), length(phenames))
Rk.coincidence=matrix(NA, length(phenames), length(phenames))
Sv.coincidence=matrix(NA, length(phenames), length(phenames))


for (l in 1:length(phenames)){ #:length(phenames)
  print(paste("-------------- Trait being analysed: ", phenames[l], " !!!!!!!!!!!---------------", sep = ""))
  
        ### get Phenotypes in rrBLUP format
        Pheno.train=as.matrix(com.tax[,phenames[l]])
        rownames(Pheno.train) <- com.tax$Taxa
        
        GEN <- Markers_impute
        GENO <- Markers_impute + 1
        M <-tcrossprod(GENO)/ncol(GENO)
        X <- GENO
        Pheno_train=Pheno.train # phenotypes
        m_train1=GEN
  
  for(j in 1:length(phenames)){


    print(paste("-------------- Training Trait: ", phenames[l], "; Validation Trait ", phenames[j], " !!!!!!!!!!!---------------", sep = ""))
    
    
    Pheno.valid=as.matrix(com.tax[,phenames[j]])
    rownames(Pheno.valid) <- com.tax$Taxa
    Pheno_valid=Pheno.valid # phenotype
    m_valid1=GEN
    
    #Step 2: Run MLMM
    Y = Pheno_train
    names(Y) <- rownames(Pheno_train)
    Y=Y[!is.na(Y)]
    
    genot.trait=m_train1[match(names(Y),rownames(m_train1)),]
    
    k2 <-Kmat[match(names(Y), rownames(Kmat)),]
    k22 <- as.matrix(k2)
    k3 <-k22[,match(names(Y), colnames(k22))]
    
    
    print(paste("-------------- Running MLMM for: ", phenames[l], "; Validation Trait ", phenames[j], " !!!!!!!!!!!---------------", sep = ""))
    
    tryCatch({
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
      
      # Fixed Effects RRBLUP
      print(paste("--------- Trait being analysed: ", phenames[l], "; Validation Trait ", phenames[j], " ; GS Model: Fixed Effects RRBLUP"," !!!!!!!----------", sep = ""))
      library(rrBLUP)
      Yt <-Pheno_train
      
      tryCatch({
        rrMod.Rx <- mixed.solve(Yt, X=m_fix_t, Z=m_train, K=NULL, SE=F, return.Hinv=F)
      }, error=function(e){})
      #tryCatch(estimatemodel(dataset), error = function() next)
      mEff.Rx<-rrMod.Rx$u
      e.Rx= as.matrix(mEff.Rx)
      
      predYv.Rx = m_valid%*%e.Rx
      
      predYr.Rx = predYv.Rx[,1]+ (m_fix_v %*% rrMod.Rx$beta)
      
      Y_valid=Pheno_valid
      Rx.accuracy[l,j] <- cor(predYr.Rx,Y_valid,use="complete")
      the.fitted.model.Rx <- lm(Y_valid ~ predYr.Rx)
      the.coefficients.Rx <- c(the.fitted.model.Rx$coefficients[1], the.fitted.model.Rx$coefficients[2])
      Rx.intercept[l,j] <- the.coefficients.Rx[1]
      Rx.slope[l,j] <- the.coefficients.Rx[2]
      
      these.observed.and.predicted.phenotypic.values.Rx <- data.frame(rownames(predYr.Rx), Y_valid, predYr.Rx)
      rownames(these.observed.and.predicted.phenotypic.values.Rx) <- NULL
      colnames(these.observed.and.predicted.phenotypic.values.Rx) <- c("Taxa", "Observed.Value", "Predicted.Value")
      x.p.Rx=these.observed.and.predicted.phenotypic.values.Rx[,c(1,3)]
      y.o.Rx=these.observed.and.predicted.phenotypic.values.Rx[,c(1,2)]
      Rx.coincidence[l,j] <- round(CI(x.p.Rx,y.o.Rx,s=s,top=T),2)
      
    }, error=function(e){})
    
    
    
    #Ordinary RRBLUP
    print(paste("--------- Trait being analysed: ", phenames[l], "; Validation Trait ", phenames[j], " ; GS Model: RRBLUP"," !!!!!!!----------", sep = ""))
    library(rrBLUP)
    Yt <-Pheno_train
    
    rrMod.Ro<-mixed.solve(Yt, X=NULL, Z=m_train1, K=NULL, SE=F, return.Hinv=F)
    mEff.Ro<-rrMod.Ro$u
    e.Ro= as.matrix(mEff.Ro)
    
    predYv.Ro = m_valid1%*%e.Ro
    
    predYr.Ro = predYv.Ro[,1]+ rrMod.Ro$beta
    
    Y_valid=Pheno_valid
    Ro.accuracy[l,j] <- cor(predYr.Ro,Y_valid,use="complete")
    the.fitted.model.Ro <- lm(Y_valid ~ predYr.Ro)
    the.coefficients.Ro <- c(the.fitted.model.Ro$coefficients[1], the.fitted.model.Ro$coefficients[2])
    Ro.intercept[l,j] <- the.coefficients.Ro[1]
    Ro.slope[l,j] <- the.coefficients.Ro[2]
    
    these.observed.and.predicted.phenotypic.values.Ro <- data.frame(names(predYr.Ro), Y_valid, predYr.Ro)
    rownames(these.observed.and.predicted.phenotypic.values.Ro) <- NULL
    colnames(these.observed.and.predicted.phenotypic.values.Ro) <- c("Taxa", "Observed.Value", "Predicted.Value")
    x.p.Ro=these.observed.and.predicted.phenotypic.values.Ro[,c(1,3)]
    y.o.Ro=these.observed.and.predicted.phenotypic.values.Ro[,c(1,2)]
    Ro.coincidence[l,j] <- round(CI(x.p.Ro,y.o.Ro,s=s,top=T),2)
    
    
    
    #RKHS 
    y.trn <- Pheno_train # for prediction accuracy
    ww <- Pheno_valid # delete data for 1/5 of the population
    
    
    print(paste("--------- Trait being analysed: ", phenames[l], "; Validation Trait ", phenames[j], " ; GS Model: RKHS"," !!!!!!!----------", sep = ""))
    ETA<-list(list(K=M,model='RKHS')) 
    fm.RK<-BGLR(y=y.trn,ETA=ETA,response_type="gaussian" ,nIter=Iter, burnIn=burn)
    Rk.accuracy[l,j] <- cor(fm.RK$yHat, ww, use="complete")
    the.fitted.model.RK <- lm(ww ~ fm.RK$yHat)
    the.coefficients.RK <- c(the.fitted.model.RK$coefficients[1], the.fitted.model.RK$coefficients[2])
    Rk.intercept[l,j] <- the.coefficients.RK[1]
    Rk.slope[l,j] <- the.coefficients.RK[2]
    
    these.observed.and.predicted.phenotypic.values.Rk <- data.frame(rownames(Y_valid), Y_valid, fm.RK$yHat)
    rownames(these.observed.and.predicted.phenotypic.values.Rk) <- NULL
    colnames(these.observed.and.predicted.phenotypic.values.Rk) <- c("Taxa", "Observed.Value", "Predicted.Value")
    x.p.Rk=these.observed.and.predicted.phenotypic.values.Rk[,c(1,3)]
    y.o.Rk=these.observed.and.predicted.phenotypic.values.Rk[,c(1,2)]
    Rk.coincidence[l,j] <- round(CI(x.p.Rk,y.o.Rk,s=s,top=T),2)
    
    
    #SVR
    print(paste("--------- Trait being analysed: ", phenames[l], "; Validation Trait ", phenames[j], " ; GS Model: SVR"," !!!!!!!----------", sep = ""))
    svp_w <- ksvm(m_train1, Yt, type="eps-svr", kernel = "rbfdot")
    yhat <- predict(svp_w, m_valid1)
    rownames(yhat) <- rownames(Y_valid)
    Sv.accuracy[l,j] <- cor(yhat, Y_valid,use="complete")
    the.fitted.model.Sv <- lm(Y_valid ~ yhat)
    the.coefficients.Sv <- c(the.fitted.model.Sv$coefficients[1], the.fitted.model.Sv$coefficients[2])
    Sv.intercept[l,j] <- the.coefficients.Sv[1]
    Sv.slope[l,j] <- the.coefficients.Sv[2]
    
    
    these.observed.and.predicted.phenotypic.values.Sv <- data.frame(rownames(Y_valid), Y_valid, yhat)
    rownames(these.observed.and.predicted.phenotypic.values.Sv) <- NULL
    colnames(these.observed.and.predicted.phenotypic.values.Sv) <- c("Taxa", "Observed.Value", "Predicted.Value")
    x.p.Sv=these.observed.and.predicted.phenotypic.values.Sv[,c(1,3)]
    y.o.Sv=these.observed.and.predicted.phenotypic.values.Sv[,c(1,2)]
    Sv.coincidence[l,j] <- round(CI(x.p.Sv,y.o.Sv,s=s,top=T),2)
    
    
    
  }

    
}


save.image("GS.Results.TraitbyTrait.BLUPs.rda")


Rx.accuracy <- data.frame(Rx.accuracy)
colnames(Rx.accuracy) <- phenames
rownames(Rx.accuracy) <- phenames
FLT.Rx <- Rx.accuracy[1:4,1:4]
write.table(Rx.accuracy, "Rx.accuracy.TraitbyTrait2.csv", sep=",", quote=F, row.names=T, col.names=T)

Ro.accuracy <- data.frame(Ro.accuracy)
colnames(Ro.accuracy) <- phenames
rownames(Ro.accuracy) <- phenames
FLT.Ro <- Ro.accuracy[1:4,1:4]
write.table(Ro.accuracy, "Ro.accuracy.TraitbyTrait2.csv", sep=",", quote=F, row.names=T, col.names=T)

Rk.accuracy <- data.frame(Rk.accuracy)
colnames(Rk.accuracy) <- phenames
rownames(Rk.accuracy) <- phenames
FLT.Rk <- Rk.accuracy[1:4,1:4]
write.table(Rk.accuracy, "Rk.accuracy.TraitbyTrait2.csv", sep=",", quote=F, row.names=T, col.names=T)


Sv.accuracy <- data.frame(Sv.accuracy)
colnames(Sv.accuracy) <- phenames
rownames(Sv.accuracy) <- phenames
FLT.Sv <- Sv.accuracy[1:4,1:4]
write.table(Sv.accuracy, "Sv.accuracy.TraitbyTrait2.csv", sep=",", quote=F, row.names=T, col.names=T)

FLT.Rx <- as.matrix(FLT.Rx)
FLT.Ro <- as.matrix(FLT.Ro)
FLT.Rk <- as.matrix(FLT.Rk)
FLT.Sv <- as.matrix(FLT.Sv)

FLT.Rx[1,1] <- 0
FLT.Rx[2,2] <- 0
FLT.Rx[3,3] <- 0
FLT.Rx[4,4] <- 0

FLT.Ro[1,1] <- 0
FLT.Ro[2,2] <- 0
FLT.Ro[3,3] <- 0
FLT.Ro[4,4] <- 0

FLT.Rk[1,1] <- 0
FLT.Rk[2,2] <- 0
FLT.Rk[3,3] <- 0
FLT.Rk[4,4] <- 0

FLT.Sv[1,1] <- 0
FLT.Sv[2,2] <- 0
FLT.Sv[3,3] <- 0
FLT.Sv[4,4] <- 0

FLT.plots <- list(FLT.Rx, FLT.Ro, FLT.Rk, FLT.Sv)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

pdf("Cowpea.FLT.Env.GS.Cor.diallele.pdf",8.2,11.2)
par(mfrow=c(4,2))
for(p in 1:4){
  par(mar=c(5,4,7,4))
    corrplot(FLT.plots[[p]], method="color", col=col(200),  
             type="upper", order="hclust", 
             addCoef.col = "black", # Add coefficient of correlation
             tl.col="black", tl.srt=45, #Text label color and rotation
             # Combine with significance
             p.mat = p.mat, sig.level = 0.05, insig = "blank", 
             # hide correlation coefficient on the principal diagonal
             diag=FALSE)
    
 
      corrplot(t(FLT.plots[[p]]), method="color", col=col(200),  
               type="upper", order="hclust", 
               addCoef.col = "black", # Add coefficient of correlation
               tl.col="black", tl.srt=45, #Text label color and rotation
               # Combine with significance
               p.mat = p.mat, sig.level = 0.05, insig = "blank", 
               # hide correlation coefficient on the principal diagonal
               diag=FALSE)
      
}
dev.off()

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#######################################################################################

library(reshape2)
models <- c("FxRRBLUP", "RRBLUP", "RKHS", "SVR")



for(p in 1:4){
  
  pdf(paste("Cowpea.FLT.Env.GS.HeatMap", models[p],"pdf", sep="."),4,4)
  
  melted_cormat <- melt(round(FLT.plots[[p]],2))
  head(melted_cormat)
  
  ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=models[p]) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1), 
          axis.text.y = element_text(vjust = 1,size = 12, hjust = 1)) + 
    coord_fixed()
  dev.off()
}

