
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
  phdata <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/Cowpea.BLUPs.csv", header = T)
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
  
  
  s <- 0.2    #  enter selection index here (proportion to select)                 
  cycles=100 # number of cross validation runs 100
  number.of.markers.to.use.for.MAS = 3
  seed.vec1 <- c(101:200)
  Rx.Trts <- NULL
  
  
  for (l in 1:length(phenames)){ #:length(phenames)
        print(paste("-------------- Trait being analysed: ", phenames[l], " !!!!!!!!!!!---------------", sep = ""))
        
        MAS.accuracy=matrix(nrow=cycles,ncol=4)
  
        
        Pheno=as.matrix(com.tax[,phenames[l]])
        rownames(Pheno) <- com.tax$Taxa
        
        for(r in 1:cycles)
        {
          print(paste("-------------- Trait being analysed: ", phenames[l], " --cycle-- ", r, " !!!!!!!!!!!---------------", sep = ""))
          #set seed set.seed(seed.vec1[r])
          set.seed(seed.vec1[r])
          #define training and validation set
          train=as.matrix(sample(1:nrow(Pheno),0.8*nrow(Pheno)), replace=F) # number of genotypes selected out of total, here 80% --> chnage 0.8 if wanted
          valid<-setdiff(1:nrow(Pheno),train)
          
          # taining set
          Pheno_train=Pheno[train,] # phenotypes
          m_train1=Markers_impute[train,]
            
          #validation set 
          Pheno_valid=Pheno[valid,] # phenotype
          m_valid1=Markers_impute[valid,]
          
          #Step 2: Run MLMM
          Y = Pheno_train
          Y=Y[!is.na(Y)]
          genot.trait=m_train1[match(names(Y),rownames(m_train1)),]
          
          k2 <-Kmat[match(names(Y), rownames(Kmat)),]
          k22 <- as.matrix(k2)
          k3 <-k22[,match(names(Y), colnames(k22))]
         
          
          print(paste("-------------- Running MLMM for: ", phenames[l], " --cycle-- ", r, " !!!!!!!!!!!---------------", sep = ""))
          
          tryCatch({
            mygwas_trait<-mlmm(Y=Y,X=genot.trait, K=k3,2,10)
          }, error=function(e){})
          #mygwas_trait<-mlmm(Y=Y,X=genot.trait, K=k3,2,10)
          res.Trait=mygwas_trait$opt_mbonf$out
          res.Trait2 <- res.Trait[order(res.Trait$pval, decreasing = F),]
          head(res.Trait2)
         # write.table(res.Trait, paste("GWAS_mlmm_results", phenames[l], "cycle",r, "csv",sep="."), sep=",", quote=F, row.names=F, col.names=T)
          
          
          #Note to Brian and Alex: we need to get this to work for multiple traits
          Gwas.output<- res.Trait2
    
          #Extract the SNP name; also extract the chromosome and bp information
          print(paste("-------Now fitting the peak marker from training set into validation set model for trait: ", phenames[l],"; cycle: ", r, " -----------", sep = ""))
          
          these.markers <- as.character(Gwas.output[1:number.of.markers.to.use.for.MAS,1])
          
          the.genotypic.data.on.these.markers <- the.genotypes[which(as.character(the.genotypes$Snp) %in% these.markers),]
          the.genotypic.data.on.these.markers.formatted.for.lm <- data.frame(colnames(the.genotypic.data.on.these.markers)[-c(1:5)],t(the.genotypic.data.on.these.markers[,-c(1:5)]))
          colnames(the.genotypic.data.on.these.markers.formatted.for.lm)[1] <- "Taxa.Names"
          colnames(the.genotypic.data.on.these.markers.formatted.for.lm)[2:4] <- these.markers
          Y.plus.order <- data.frame(Pheno,1:nrow(Pheno))
          Y.plus.order$Taxa.Names <- rownames(Y.plus.order)
          colnames(Y.plus.order)[1] <- "Trait"
          
          data.for.lm.almost <- merge(Y.plus.order, the.genotypic.data.on.these.markers.formatted.for.lm, by = "Taxa.Names")
          data.for.lm.almost <- data.for.lm.almost[order(data.for.lm.almost[,3]),]
        
          
       
          data.for.lm.train <- data.for.lm.almost[train,]
          data.for.lm.pred <- data.for.lm.almost[valid,]
          equation.for.lm <- paste(colnames(data.for.lm.train)[2], "~", colnames(data.for.lm.train)[4],sep = "")
          for(index in 5:ncol(data.for.lm.train)) equation.for.lm <- paste(equation.for.lm,colnames(data.for.lm.train)[index],sep = "+")
          
          
          lm.model.fitted.in.train <- lm(equation.for.lm, data = data.for.lm.train)
          the.predicted.MAS.values <- predict(lm.model.fitted.in.train, newdata = data.for.lm.pred)
          
          these.observed.and.predicted.phenotypic.values <- data.frame(data.for.lm.pred[,1:2], the.predicted.MAS.values)
          colnames(these.observed.and.predicted.phenotypic.values) <- c("Taxa", "Observed.Value", "Predicted.Value")
          
          #Measure correclation between OBS and Pred in validation set (V.S.)
          MAS.accuracy[r,1] <- cor(these.observed.and.predicted.phenotypic.values[,3], these.observed.and.predicted.phenotypic.values[,2],use="complete")
          x.p.MAS=these.observed.and.predicted.phenotypic.values[,c(1,3)]
          y.o.MAS=these.observed.and.predicted.phenotypic.values[,c(1,2)]
          MAS.accuracy[r,4] <- round(CI(x.p.MAS,y.o.MAS,s=s,top=T),2)
          #Fit a linear regression model where the observed values are the response variable and the predicted values is the explanatory variable 
            the.fitted.regression.model <- lm(these.observed.and.predicted.phenotypic.values[,2] ~ these.observed.and.predicted.phenotypic.values[,3])
          
          
          
          MAS.accuracy[r,2] <- the.fitted.regression.model$coefficients[1]
          MAS.accuracy[r,3] <- the.fitted.regression.model$coefficients[2]
  
        }
          print(paste("--------- Summarizing Result for trait: ", phenames[l]," !!!!!!!----------", sep = ""))
          MAS.accuracy <- data.frame(MAS.accuracy)
          MAS.accuracy$Trt <- rep(phenames[l], nrow(MAS.accuracy))
          Rx.Trts <- rbind(Rx.Trts, MAS.accuracy)
        
    }

  save.image("MAS.Results.100cycles.BLUPs.rda")
  Rx.Trts$Model <- rep("MAS", nrow(Rx.Trts))
  colnames(Rx.Trts)[c(1:4)] <- c("P.Accuracy", "Intercept", "Slope", "C.Index")

  MAS.Results <- Rx.Trts
  write.table(MAS.Results, "Vigna.MAS.Results.AllTraits.BLUPs.12232018.txt", sep="\t", quote=F, row.names = F, col.names = T)  

  #colnames(Rx.Trts)[1:4] <- c("X1", "X2", "X3", "X4")
  
    