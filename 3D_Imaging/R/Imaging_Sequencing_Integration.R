library(matlab)
library(ggplot2)
setwd("~/Projects/PMO/MeasuringFitnessPerClone/code/SingleCellSequencing")
ROOT="~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87"
OUTCORRECTED=paste0(ROOT,filesep,"A05_PostProcessCellposeOutput")
K=15 ## neighbors
N=30; ## cells
MINNUCVOL=8^3
FoF="FoF16_210818_fluorescent.cytoplasm"
# FoF="FoF13_220228_fluorescent.cytoplasm"

overlayHist<-function(this,ontothat){
  dat=as.data.frame(c(this,ontothat))
  colnames(dat)="data"
  dat$cat="ontothat"
  dat$cat[1:length(this)]="this"
  p=ggplot(dat,aes(x=data,fill=cat)) + 
    geom_histogram(data=dat, alpha = 0.4)
  return(list(dat=dat,p=p))
}
# this=imgStats[,1]; ontothat=sample(imgStats[,2],100)
overlayDistributions1D<-function(this, ontothat,q=c(0.05,0.5,0.95)){
  # dat=overlayHist(this,ontothat)
  qstat<-function(this, ontothat){
    q=sapply(list(this,ontothat), quantile, q)
    colnames(q)=c("this","ontothat")
    q=as.data.frame(q)
    return(q)
  }
  q1=qstat(this, ontothat)
  ## Shift both medians to 0
  this=this-q1$this[1]
  ontothat=ontothat-q1$ontothat[1]
  ## Overlay upper quantile
  q2=qstat(this, ontothat)
  this=this/(q2$this[3]/q2$ontothat[3])
  
  ## ontothat: Median back to orig
  ontothat=ontothat+q1$ontothat[2]
  this=this+q1$ontothat[2]
  dat=overlayHist(this,ontothat)
  return(dat)
}



## read seq stats
# seqStats=read.table("../../data/GastricCancerCL/RNAsequencing/B02_220112_seqStats/NCI-N87/Clone_0.244347_ID119967.txt",sep="\t",check.names = F,stringsAsFactors = F)
load('~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/RNAsequencing/B01_220112_pathwayActivity/NCI-N87/Clone_0.244347_ID119967.RObj')
ccState=sapply(colnames(pq), function(x) cloneid::getAttribute(x,"TranscriptomePerspective","state"))
seqStats=t(pq)

## read imaging stats
imgStats=read.table(paste0("../../data/GastricCancerCL/3Dbrightfield/NCI-N87/A07_LinkedSignals_Stats/",FoF,"_stats.txt"),sep="\t",check.names = F,stringsAsFactors = F,header = T)
imgStats=imgStats[apply(!is.na(imgStats),1,all),]; ## exclude cells whose volume could not be estimated
imgStats=imgStats[,apply(!is.na(imgStats),2,all)]; ## exclude features with NA vals
imgStats=imgStats[imgStats$vol_nucleus>MINNUCVOL,]

###################################
## Compare imaging and seq stats ##
###################################
## Visualize spatial distribution of cell IDs
centers=read.csv(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates/nucleus.p_Cells_Centers.csv"))
centers=centers[as.character(rownames(imgStats)),]
## Hypothesis: cells in the center of islands don't divide -- top k dense cells are considered G0/G1 cells
d=dist(centers)
d=apply(as.matrix(d),1, function(x) mean(sort(x)[1:K]))
g1cells_img=names(sort(d)[1:N])
plot(centers$x,-centers$y,cex=log(imgStats[,"area_nucleus.p"]/1000),col=1+rownames(centers) %in% g1cells_img,pch=1+19*rownames(centers) %in% g1cells_img)
text(centers$x,-centers$y, rownames(centers),cex=0.4)

## Image stats postprocess
# imgStats=as.data.frame(apply(imgStats, 2, function(x) (x - mean(x,na.rm=T))/sd(x,na.rm=T))) ## standardize
# write.table(cbind(rownames(imgStats),imgStats),col.names = c("cellID",colnames(imgStats)),file=paste0("~/Downloads/",FoF,"_organelleFeatures.txt"),sep="\t",quote = F,row.names = F)

## Same number of sequenced and imaged cells:
seqStats=seqStats[sample(nrow(seqStats),size = nrow(imgStats)),]
g1cells_seq=rownames(seqStats)[ccState[rownames(seqStats)]=="G0G1"]

## UMAP
imgStats_=as.data.frame(umap::umap(imgStats)$layout)
seqStats_=as.data.frame(umap::umap(seqStats)$layout)

## MST imaging:
par(mfrow=c(1,2))
tree=slingshot(imgStats_,rep(1,nrow(imgStats_))); 
cc_img=slingshot::slingPseudotime(tree)+1
col=rainbow(max(cc_img[rownames(imgStats_),1]))
plot(imgStats_, asp = 1,pch=15-13*(!rownames(imgStats_) %in% g1cells_img), col=col[round(cc_img[rownames(imgStats_),1])])
lines(as.SlingshotDataSet(tree), type = 'c', lwd = 3)
## MST sequencing:
tree2=slingshot(seqStats_,rep(1,nrow(seqStats_))); 
cc_seq=slingshot::slingPseudotime(tree2)+1
col=rainbow(max(cc_seq[rownames(seqStats_),1]))
plot(seqStats_, asp = 1, pch=15-13*(!rownames(seqStats_) %in% g1cells_seq),col=col[round(cc_seq[rownames(seqStats_),1])])
lines(as.SlingshotDataSet(tree2), type = 'c', lwd = 3)

## Calculate pseudotime difference to G1 cells
d_i= sapply(cc_img, function(x) quantile(cc_img[g1cells_img,]-x,))
d_s= sapply(cc_seq, function(x) quantile(cc_seq[g1cells_seq,]-x,))
d_i= as.data.frame(t(d_i))
d_s= as.data.frame(t(d_s))
d_i$type="img"
d_s$type="seq"

## co-cluster image and sequencing stats
stats=rbind(d_i,d_s)
rownames(stats) = paste(rownames(stats), stats$type)
dd = dist(stats[,1:2])
tr = ape::nj(dd)
col = rep("red", length(tr$tip.label))
col[grep("img",tr$tip.label) ] = "blue"
par(mfrow=c(1,1))
plot(tr,show.tip.label = T, tip.color = col, cex=0.36)
legend("topright",c("sequenced cell","imaged cell"),fill=c("red","blue"),cex=1.8)
