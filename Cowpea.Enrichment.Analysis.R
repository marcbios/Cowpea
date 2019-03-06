rm(list=ls())

setwd("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Genes/")

ft.genes <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Genes/Cowpea.AprioriGenes.02112019.csv", header=T)
head(ft.genes, 30)
epi.QTL <- read.csv("Cowpea.EpiQTL.PhysicalPosition.csv", header=T)
tail(ft.genes)
jl.QTL <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.JL.QTL.PhysicalPos.csv", header=T)

str(ft.genes)
ft.genes$Gene.name <- as.character(ft.genes$Gene.name)
ft.genes$Flowering.pathway <- as.character(ft.genes$Floweringpathway)
ft.genes$Cowpea.GeneID <- as.character(ft.genes$Cowpea.GeneID)

str(jl.QTL)
jl.QTL$QTL <- as.character(jl.QTL$QTL)
jl.QTL$Trait <- as.character(jl.QTL$Trait)
jl.QTL$QTL2 <- as.character(jl.QTL$QTL2)

write.table(t(data.frame(c("QTL", "LG", "Position", "Trait", "PVE", "QTL2", names(ft.genes), "gene.pres", "Q.Dist1", "Q.Dist2", "M.Dist1", "M.Dist2"))),paste("Cowpea.QTL_GeneEnrichment.02172019", ".csv", sep=""), sep=",", append=T, quote=F, row.names=F, col.names=F)

for(q in 1:nrow(jl.QTL)){
    
    for(g in 1:nrow(ft.genes)){
      
        if(jl.QTL[q,9]==ft.genes[g,8]){
          
            if(ft.genes[g,9]>=jl.QTL[q,10]&ft.genes[g,9]<=jl.QTL[q,11]){
              gene.pres <- "Y"
            }else(gene.pres <- "N")
            
            POS1 <- jl.QTL[q,10] - ft.genes[g,9]
            POS2 <- jl.QTL[q,11] - ft.genes[g,10]
            POS3 <- jl.QTL[q,13] - ft.genes[g,9]
            POS4 <- jl.QTL[q,14] - ft.genes[g,9]
            
            gene_details <- as.vector(as.matrix(ft.genes[g, c(1:10)]))
            
            write.table(t(data.frame(c(jl.QTL[q,1], jl.QTL[q,2], jl.QTL[q,3], jl.QTL[q,4], jl.QTL[q,6], jl.QTL[q,8], gene_details, gene.pres[1], POS1[1], POS2[1], POS3[1], POS4[1]))),paste("Cowpea.QTL_GeneEnrichment.02172019", ".csv", sep=""), sep=",", append=T, quote=F, row.names=F, col.names=F)
            
        }else(next)
      
    }
  
  
}



## Epistasis QTL
epi.qtl <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/SPAEML.Epistatic.QTL.Pos.Summary.New.02.09.2019.csv", header = T)

head(epi.qtl)
epi.qtl.1 <- epi.qtl[,c(1,2,4:9)]
epi.qtl.2 <- epi.qtl[,c(4,10:16)]

names(epi.qtl.1)
names(epi.qtl.2)

names(epi.qtl.1) <- c("Marker", "QTL", "Trait", "LG", "Position", "Additive.Effect", "PVE","MAF")

names(epi.qtl.2) <- c("Trait", "Marker", "QTL", "LG", "Position","Additive.Effect", "PVE", "MAF")

epi.qtl.1$Epi <- rep("LOC1", nrow(epi.qtl.1))
epi.qtl.2$Epi <- rep("LOC2", nrow(epi.qtl.2))

epi.qtl.both <- rbind(epi.qtl.1, epi.qtl.2)
head(epi.qtl.both)

jl.QTL <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.JL.QTL.PhysicalPos.csv", header=T)
head(jl.QTL)
names(jl.QTL)[1] <- "Marker"
names(jl.QTL)[8] <- "QTL"
sel.qtl <- jl.QTL[,c(1,8,9:14)]
head(sel.qtl)
epi.qtl.both2 <- merge(epi.qtl.both, sel.qtl, all.x = T)
write.table(epi.qtl.both2, "/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Cowpea.EpiQTL.PhysicalPos.csv", sep=",", quote=F, row.names = F, col.names = T)



rm(list=ls())

setwd("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Genes/")

ft.genes <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Genes/Cowpea.AprioriGenes.02112019.csv", header=T)
head(ft.genes, 30)

tail(ft.genes)
epi.QTL <- read.csv("/Users/omo/Google Drive/Post Doc/Collaborative Publications/Cowpea/Genes/Cowpea.EpiQTL.PhysicalPosition.csv", header=T)

str(ft.genes)
ft.genes$Gene.name <- as.character(ft.genes$Gene.name)
ft.genes$Flowering.pathway <- as.character(ft.genes$Floweringpathway)
ft.genes$Cowpea.GeneID <- as.character(ft.genes$Cowpea.GeneID)

str(epi.QTL)
epi.QTL$QTL <- as.character(epi.QTL$QTL)
epi.QTL$Trait <- as.character(epi.QTL$Trait)
epi.QTL$QTL <- as.character(epi.QTL$QTL)
epi.QTL$QTL <- as.character(epi.QTL$QTL)

write.table(t(data.frame(c("Marker","QTL", "LG", "Position", "Trait", "LOC","PVE", names(ft.genes), "gene.pres", "Q.Dist1", "Q.Dist2"))),paste("Cowpea.Epi.QTL_GeneEnrichment.02182019", ".csv", sep=""), sep=",", append=T, quote=F, row.names=F, col.names=F)

for(q in 1:nrow(epi.QTL)){
  
  for(g in 1:nrow(ft.genes)){
    
    if(epi.QTL[q,10]==ft.genes[g,8]){
      
      if(ft.genes[g,9]>=epi.QTL[q,11]&ft.genes[g,9]<=epi.QTL[q,12]){
        gene.pres <- "Y"
      }else(gene.pres <- "N")
      
      POS1 <- epi.QTL[q,11] - ft.genes[g,9]
      POS2 <- epi.QTL[q,12] - ft.genes[g,10]
      
      gene_details <- as.vector(as.matrix(ft.genes[g, c(1:10)]))
      
      write.table(t(data.frame(c(as.character(epi.QTL[q,1]), epi.QTL[q,2], epi.QTL[q,4], epi.QTL[q,5], epi.QTL[q,3], epi.QTL[q,9], epi.QTL[q,7], gene_details, gene.pres[1], POS1[1], POS2[1]))),paste("Cowpea.Epi.QTL_GeneEnrichment.02182019", ".csv", sep=""), sep=",", append=T, quote=F, row.names=F, col.names=F)
      
    }else(next)
    
  }
  
  
}


