library(matlab)
library(flexclust)
library(ggplot2)
library(slingshot)
setwd("~/Projects/PMO/MeasuringFitnessPerClone/code/SingleCellSequencing")
ROOT="~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87"
INCORRECTED=paste0(ROOT,filesep,"A05_PostProcessCellposeOutput")
INSTATS=paste0(ROOT,filesep,"A07_LinkedSignals_Stats")
OUTPSEUDOTIME=paste0(ROOT,filesep,"A08_Pseudotime")
dirCreate(OUTPSEUDOTIME,permission = "a+w")
FROMILASTIK=paste0(ROOT,filesep,"G07_IlastikOutput")
xyz=c("x","y","z")
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

## read imaging stats (if timeseries, FoFs must be sorted in ascending temporal order)
# FoFs=c("FoF13_220228_fluorescent.cytoplasm","FoF16_210818_fluorescent.cytoplasm")
FoFs=paste0("FoF",1:4,"007_220523_brightfield")
imgStats=list()
for(FoF in FoFs){
  imgStats_=read.table(paste0(INSTATS,filesep,FoF,"_stats.txt"),sep="\t",check.names = F,stringsAsFactors = F,header = T)
  imgStats_=imgStats_[apply(!is.na(imgStats_),1,all),]; ## exclude cells whose volume could not be estimated
  imgStats_=imgStats_[,apply(!is.na(imgStats_),2,all)]; ## exclude features with NA vals
  imgStats_$FoF= FoF
  ## which frame (timepoint) is this? Used for mapping to Ilastik output:
  imgStats_$frame=which(FoF==FoFs)-1;  
  imgStats[[FoF]]=imgStats_[imgStats_$vol_nucleus>MINNUCVOL,]
}
# run pseudotime inference on all timepoints of a given FoV jointly:
imgStats=do.call(rbind, imgStats)
print(paste("Found",nrow(imgStats),"cells across",length(FoFs),"images"))
## Columns of interest for slingshot
coi=setdiff(colnames(imgStats),c("FoF","frame")); #xyz

###################################
## Compare imaging and seq stats ##
###################################
## Visualize spatial distribution of cell IDs
centers=read.csv(paste0(INCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates/nucleus.p_Cells_Centers.csv"))
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
imgStats_=as.data.frame(umap::umap(apply(imgStats[,coi],2,as.numeric))$layout)
seqStats_=as.data.frame(umap::umap(seqStats)$layout)

## MST imaging:
pdf("~/Downloads/slingshot_imaging.pdf")
par(mfrow=c(2,2))
tree=slingshot(imgStats_,rep(1,nrow(imgStats_))); 
pseudotime_img=as.data.frame(slingshot::slingPseudotime(tree)+1)
pseudotime_img$col=round(pseudotime_img[rownames(imgStats_),1])
col=imgStats$frame+1
plot(imgStats_, asp = 1,pch=15-13, col=col)
lines(as.SlingshotDataSet(tree), type = 'c', lwd = 3)
color.bar(unique(col),min=1,max = max(col),nticks = length(unique(col)),title = "real order")
col=c(brewer.pal(9,"Reds")[3:9],brewer.pal(9,"YlGn")[3:9],brewer.pal(9,"Blues")[3:9])[1:max(pseudotime_img$col)]
# plot(imgStats_, asp = 1,pch=15-13*(!rownames(imgStats_) %in% g1cells_img), col=col[pseudotime_img$col])
plot(imgStats_, asp = 1,pch=2, col=col[pseudotime_img$col])
lines(as.SlingshotDataSet(tree), type = 'c', lwd = 3)
color.bar(col,min=1,max = length(col),title = "pseudo order")
dev.off()
## MST sequencing:
tree2=slingshot(seqStats_,rep(1,nrow(seqStats_))); 
pseudotime_seq=slingshot::slingPseudotime(tree2)+1
plot(seqStats_, asp = 1, pch=15-13*(!rownames(seqStats_) %in% g1cells_seq),col=col[round(pseudotime_seq[rownames(seqStats_),1])])
lines(as.SlingshotDataSet(tree2), type = 'c', lwd = 3)


## Calculate pseudotime difference to G1 cells
d_i= sapply(pseudotime_img, function(x) quantile(pseudotime_img[g1cells_img,]-x,))
d_s= sapply(pseudotime_seq, function(x) quantile(pseudotime_seq[g1cells_seq,]-x,))
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


## Save for comparison with live-cell imaging data
pseudotime_img=as.data.frame(pseudotime_img)
pseudotime_img[,c(xyz,"FoF","frame")]=imgStats[,c(xyz,"FoF","frame")]
FoV=gsub(substr(FoF,1,4), "FoFX",FoF)
write.table(pseudotime_img, file=paste0(OUTPSEUDOTIME,filesep,FoV,".txt"),sep="\t",row.names = F,quote = F)



#############################################################
## Compare pseudotime with LCI-derived time since division ##
#############################################################
realtime_img=read.table(file=paste0(FROMILASTIK,filesep,FoF,"_DeltaDivision.txt"),sep="\t",check.names = T,stringsAsFactors = T,header = T)
## Coordinate based mapping of cells between pseudotime and Ilastik output per each frame:
jointTimes=list()
par(mfrow=c(3,2))
for(time in unique(pseudotime_img$frame)){
  print(paste("merging time",time,"..."))
  pseudotime_img_ = pseudotime_img[pseudotime_img$frame==time,]
  realtime_img_ = realtime_img[realtime_img$frame==time,]
  d=dist2(pseudotime_img_[,xyz], realtime_img_[,xyz])
  ii=apply(d,2,which.min)
  realtime_img_$pseudotime= pseudotime_img_$Lineage1[ii] 
  ##@TODO: test -- there should be a unique arg minimum for each cell
  print(apply(pseudotime_img[,xyz],2,quantile))
  print(apply(realtime_img_[,xyz],2,quantile))
  # heatmap.2(d,trace="n")
  hist(plyr::count(ii)$freq,30,main=paste("time",time))
  # - focus on just cells with assigned parent track ID
  ii=!is.na(realtime_img_$time_since_division) | !is.na(realtime_img_$time_until_division)
  jointTimes[[as.character(time)]]=realtime_img_[ii,,drop=F]
}
jointTimes=do.call(rbind, jointTimes)
write.csv(jointTimes,file=paste0(OUTPSEUDOTIME,filesep,FoV,".csv"),row.names = F)
# # run in Matlab: 
#   # calc circular cross correlation btw. the "pseudotime timepoints" and the estimates of "time_since_division"
#   circularCrossCorr(jointTimes$time_since_division,jointTimes$pseudotime)
## Read in data after shifting pseudotime by cross correlation
jointTimes=read.csv(file=paste0(OUTPSEUDOTIME,filesep,FoV,"_matlabOut.csv"),check.names = F,stringsAsFactors = F)
te=cor.test(jointTimes$time_since_division,jointTimes$pseudotime_shifted,method = "spearman")
pdf("~/Downloads/realtime_vs_pseudotime.pdf",width = 4.5,height = 4.5)
plot(jointTimes$time_since_division,jointTimes$pseudotime_shifted,col=col[jointTimes$pseudotime],pch=20,xlab="real time since division (hours)",ylab="pseudotime",cex=2,main=paste0("r=",round(te$estimate,2),"; P=",round(te$p.value,2)))
# points(jointTimes$time_since_division,jointTimes$pseudotime_shifted,cex=2)
dev.off()


