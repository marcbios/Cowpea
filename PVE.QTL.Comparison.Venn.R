QTL.EPI <- read.csv("~/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Data/TASSEL/NewData/Main.EPI.QTL.csv", header=T)
head(QTL.EPI)
QTL.EPI.FT <- QTL.EPI[which(QTL.EPI$Trait=="FLT"),]
head(QTL.EPI.FT)


source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")

setlist <- list(Q.MAP=as.vector(QTL.EPI.FT$QTL1[which(QTL.EPI.FT$Epi=="Q.MAP")]), E.MAP=as.vector(QTL.EPI.FT$QTL1[which(QTL.EPI.FT$Epi=="E.MAP")]))

setlist2 <- setlist[1:2]
OLlist2 <- overLapper(setlist=setlist2, sep="_", type="vennsets")
counts <- list(sapply(OLlist2$Venn_List, length))
pdf("FLT_Venn.pdf")
vennPlot(counts=counts, mysub="Flowering Time", yoffset=c(0.3, -0.2)) 
dev.off()


QTL.EPI.MT <- QTL.EPI[which(QTL.EPI$Trait=="MAT"),]
head(QTL.EPI.MT)

setlist <- list(Q.MAP=as.vector(QTL.EPI.MT$QTL1[which(QTL.EPI.MT$Epi=="Q.MAP")]), E.MAP=as.vector(QTL.EPI.MT$QTL1[which(QTL.EPI.MT$Epi=="E.MAP")]))

setlist2 <- setlist[1:2]
OLlist2 <- overLapper(setlist=setlist2, sep="_", type="vennsets")
counts <- list(sapply(OLlist2$Venn_List, length))
pdf("MAT_Venn.pdf")
vennPlot(counts=counts, mysub="Flowering Time", yoffset=c(0.3, -0.2)) 
dev.off()



QTL.EPI.SS <- QTL.EPI[which(QTL.EPI$Trait=="SSZ"),]
head(QTL.EPI.SS)

setlist <- list(Q.MAP=as.vector(QTL.EPI.SS$QTL1[which(QTL.EPI.SS$Epi=="Q.MAP")]), E.MAP=as.vector(QTL.EPI.SS$QTL1[which(QTL.EPI.SS$Epi=="E.MAP")]))

setlist2 <- setlist[1:2]
OLlist2 <- overLapper(setlist=setlist2, sep="_", type="vennsets")
counts <- list(sapply(OLlist2$Venn_List, length))

pdf("SS_Venn.pdf")
vennPlot(counts=counts, mysub="Flowering Time", yoffset=c(0.3, -0.2)) 
dev.off()













