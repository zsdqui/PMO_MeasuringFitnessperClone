source("~/Projects/code/RCode/scripts/grpstats.R")

setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/LiveCellImaging/A03_190516_Tracking_ilastik")
wells_96=matrix(NA, 8, 12)
rownames(wells_96)= c("A","B","C","D","E","F","G","H")
colnames(wells_96)=c(paste0("0",1:9),10:12)
MAXCELLSPERWELL = 120
MINCELLSPERWELL = 60

wells_96 = wells_96[1:3,1:3]; ##Subset

pdf(file = "~/Downloads/countePerWell.pdf", width = 6, height = 6)
par(mgp = c(2, 1, 0),mfrow=c(nrow(wells_96), ncol(wells_96)), mai = c(0.4, 0.4, 0.3, 0.1))
for(i in rownames(wells_96)){
  for(j in colnames(wells_96)){
    dm=read.table(paste0(i,j,"_CellTracks_CSV-Table.csv"),sep=",", header = T, stringsAsFactors = F, check.names = F)
    dm=dm[dm$trackId!=-1,]
    cellcnt=sapply(unique(dm$frame), function(x) length(unique(dm$trackId[dm$frame==x])))
    names(cellcnt)=paste0(1:12,"h")
    wells_96[i,j] = cellcnt[length(cellcnt)]/cellcnt[1]
    plot(cellcnt, xlab="hour", ylab="cell count", ylim=c(MINCELLSPERWELL,MAXCELLSPERWELL), log="y", pch = 20, cex=2)
    title(paste("GR =",round(wells_96[i,j],2)), line=0.5)
  }
}
dev.off()


