library(matlab)
library(RColorBrewer)
library(flexclust)
library(ggplot2)
library(slingshot)
source("~/Projects/code/RCode/scripts/color.bar.R")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
setwd("~/Projects/PMO/MeasuringFitnessPerClone/code/SingleCellSequencing")
ROOT="~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87"
A01=paste0(ROOT,filesep,"A01_rawData")
INCORRECTED=paste0(ROOT,filesep,"A05_PostProcessCellposeOutput")
INSTATS=paste0(ROOT,filesep,"A07_LinkedSignals_Stats")
OUTPSEUDOTIME=paste0(ROOT,filesep,"A08_Pseudotime")
dirCreate(OUTPSEUDOTIME,permission = "a+w")
FROMILASTIK=paste0(ROOT,filesep,"G07_IlastikOutput")
xyz=c("x","y","z")
K=15 ## neighbors
N=30; ## cells
MINNUCVOL=8^3
xmlfiles=list.files('../../data/GastricCancerCL/3Dbrightfield/NCI-N87/A01_rawData/',pattern=".xml",full.names=T)


## read seq stats
# seqStats=read.table("../../data/GastricCancerCL/RNAsequencing/B02_220112_seqStats/NCI-N87/Clone_0.244347_ID119967.txt",sep="\t",check.names = F,stringsAsFactors = F)
load('~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/RNAsequencing/B01_220112_pathwayActivity/NCI-N87/Clone_0.244347_ID119967.RObj')
ccState=sapply(colnames(pq), function(x) cloneid::getAttribute(x,"TranscriptomePerspective","state"))
seqStats=t(pq)

## read imaging stats (if timeseries, FoFs must be sorted in ascending temporal order)
# FoFs=grep("00200", gsub("_stats.txt","",list.files(INSTATS,pattern = "_221018_brightfield")), value=T)
FoFs=grep("00100", gsub("_stats.txt","",list.files(INSTATS,pattern = "_221018_brightfield")), value=T)
# FoFs=paste0("FoF",1:5,"003_220721_brightfield")
imgStats=list()
for(FoF in FoFs){
  imgStats_=read.table(paste0(INSTATS,filesep,FoF,"_stats.txt"),sep="\t",check.names = F,stringsAsFactors = F,header = T)
  imgStats_=imgStats_[apply(!is.na(imgStats_),1,all),]; ## exclude cells whose volume could not be estimated
  imgStats_=imgStats_[,apply(!is.na(imgStats_),2,all)]; ## exclude features with NA vals
  imgStats_$FoF= FoF
  ## which frame (timepoint) is this? Used for mapping to Ilastik output:
  imgStats_$frame=which(FoF==FoFs)-1;  
  imgStats_$ID=rownames(imgStats_)
  imgStats[[FoF]]=imgStats_[imgStats_$vol_nucleus>MINNUCVOL,]
}
# run pseudotime inference on all timepoints of a given FoV jointly:
imgStats=do.call(rbind, imgStats)
print(paste("Found",nrow(imgStats),"cells across",length(FoFs),"images"))
## Columns of interest for slingshot
coi=setdiff(colnames(imgStats),c("FoF","frame","ID","count_nucleus.p",xyz)); #
# coi=grep("mito",coi,value=T,invert = T)
imgStats[,coi]=sweep(imgStats[,coi], 2, STATS = apply(imgStats[,coi],2,median),FUN = "/")

## Same number of sequenced and imaged cells:
# seqStats=seqStats[sample(nrow(seqStats),size = nrow(imgStats)),]
# g1cells_seq=rownames(seqStats)[ccState[rownames(seqStats)]=="G0G1"]

## UMAP
imgStats_=as.data.frame(umap::umap(apply(imgStats[,coi],2,as.numeric))$layout)
# seqStats_=as.data.frame(umap::umap(seqStats)$layout)

## MST imaging:
tree=slingshot(imgStats_,rep(1,nrow(imgStats_))); 
pseudotime_img=as.data.frame(slingshot::slingPseudotime(tree)+1)
colnames( pseudotime_img)="pseudotime"
## Save for comparison with live-cell imaging data
pseudotime_img=as.data.frame(pseudotime_img)
pseudotime_img[,c(xyz,"FoF","frame","cellID")]=imgStats[,c(xyz,"FoF","frame","ID")]
hours=getTimeStampsFromMetaData(FoFs, root="A01_rawData", xmlfiles)
pseudotime_img$hour=hours[match(pseudotime_img$FoF,names(hours))]
FoV=gsub(substr(FoFs[1],1,4), "FoFX",FoFs[1])
# write.table(pseudotime_img, file=paste0(OUTPSEUDOTIME,filesep,FoV,".txt"),sep="\t",row.names = F,quote = F)
write.csv(pseudotime_img,   file=paste0(OUTPSEUDOTIME,filesep,FoV,".csv"),row.names = F)


## Visualize pseudotime
pdf(paste0("~/Downloads/",FoV,"_slingshot_imaging.pdf"))
par(mfrow=c(2,2))
pseudotime_img$col=round(pseudotime_img[rownames(imgStats_),1])
col=rainbow(max(pseudotime_img$hour)*1.2)
col=col[round(pseudotime_img$hour)]
plot(imgStats_, asp = 1,pch=20, col=col)
lines(as.SlingshotDataSet(tree), type = 'c', lwd = 3)
color.bar(unique(col),min=min(round(pseudotime_img$hour)),max = max(round(pseudotime_img$hour)),nticks = length(unique(pseudotime_img$hour)),title = "real order")
col=c(brewer.pal(9,"Reds")[3:9],brewer.pal(9,"YlGn")[3:9],brewer.pal(9,"Blues")[3:9])[1:max(pseudotime_img$col)]
# plot(imgStats_, asp = 1,pch=15-13*(!rownames(imgStats_) %in% g1cells_img), col=col[pseudotime_img$col])
plot(imgStats_, asp = 1,pch=20, col=col[pseudotime_img$col])
lines(as.SlingshotDataSet(tree), type = 'c', lwd = 3)
color.bar(col,min=1,max = length(col),title = "pseudo order")
dev.off()
print("Correlation between real time order and pseudotime:")
cor.test(pseudotime_img$frame,pseudotime_img$pseudotime,method="spearman")
print("Correlation between real time and pseudotime:")
te=cor.test(pseudotime_img$hour,pseudotime_img$pseudotime)
boxplot(pseudotime_img$pseudotime~round(pseudotime_img$hour),main=paste0("Pearson r=",round(te$estimate,2),"; P=",round(te$p.value,5)),col="cyan",xlab="hour",ylab="pseudotime")



## Visualize pseudotime spatial distribution
par(mfrow=c(5,3),bg="gray",mai=c(0.1,0.5,0.1,0.5))
zslice=35
pseudotime_img[,"pseudotimeCol"]=round(pseudotime_img[,"pseudotime"]*10)
pseudotime_img[,"randCol"]=sample(pseudotime_img[,"pseudotimeCol"])
col=fliplr(heat.colors(max(pseudotime_img[,"pseudotimeCol"])))
for(FoF in unique(pseudotime_img$FoF)){
  pseudotime_img_=pseudotime_img[pseudotime_img$FoF==FoF,]
  rownames(pseudotime_img_)=as.character(pseudotime_img_$cellID)
  slice=getZslice(FoF,zslice,root = "../../data/GastricCancerCL/3Dbrightfield/NCI-N87/A06_multiSignals_Linked",plot=F)
  ## plot brightfield
  TIF=paste0("../../data/GastricCancerCL/3Dbrightfield/NCI-N87/A01_rawData",filesep,FoF,filesep,"nucleus.s_z",zslice,".tif")
  img=bioimagetools::readTIF(TIF)
  img=EBImage::rotate(img,-180)
  bioimagetools::img(resize4Ilastik(img, xydim = xydim)[,,1]);
  ## plot pseudotime colorcoded
  plot(-slice$x,slice$y,col=col[pseudotime_img_[as.character(slice$ID),"pseudotimeCol"]],main=FoF,xaxt='n',yaxt='n', ann=FALSE)
  ## plot random color code
  plot(-slice$x,slice$y,col=col[pseudotime_img_[as.character(slice$ID),"randCol"]],main=FoF,xaxt='n',yaxt='n', ann=FALSE)
}
color.bar(unique(col),min=1,max = length(col),nticks = length(unique(col)),title = "pseudo order")
## zoom in:
ii=which(slice$x>0.5*max(slice$x) )
plot(-slice$x[ii],slice$y[ii],col=col[pseudotime_img_[as.character(slice$ID[ii]),"pseudotimeCol"]],main=FoF,xaxt='n',yaxt='n', ann=FALSE,ylim=quantile(slice$y,c(0,1))*c(1,0.75))
file.copy(TIF,"~/Downloads/")
print(paste("open",TIF))


# ## MST sequencing:
# tree2=slingshot(seqStats_,rep(1,nrow(seqStats_))); 
# pseudotime_seq=slingshot::slingPseudotime(tree2)+1
# plot(seqStats_, asp = 1, pch=15-13*(!rownames(seqStats_) %in% g1cells_seq),col=col[round(pseudotime_seq[rownames(seqStats_),1])])
# lines(as.SlingshotDataSet(tree2), type = 'c', lwd = 3)
# 
# 
# ## Calculate pseudotime difference to G1 cells
# d_i= sapply(pseudotime_img, function(x) quantile(pseudotime_img[g1cells_img,]-x,))
# d_s= sapply(pseudotime_seq, function(x) quantile(pseudotime_seq[g1cells_seq,]-x,))
# d_i= as.data.frame(t(d_i))
# d_s= as.data.frame(t(d_s))
# d_i$type="img"
# d_s$type="seq"
# 
# ## co-cluster image and sequencing stats
# stats=rbind(d_i,d_s)
# rownames(stats) = paste(rownames(stats), stats$type)
# dd = dist(stats[,1:2])
# tr = ape::nj(dd)
# col = rep("red", length(tr$tip.label))
# col[grep("img",tr$tip.label) ] = "blue"
# par(mfrow=c(1,1))
# plot(tr,show.tip.label = T, tip.color = col, cex=0.36)
# legend("topright",c("sequenced cell","imaged cell"),fill=c("red","blue"),cex=1.8)
# 


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
  realtime_img_$pseudotime= pseudotime_img_$pseudotime[ii] 
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
# te=cor.test(jointTimes$time_since_division,jointTimes$pseudotime_shifted,method = "spearman")
te=cor.test(jointTimes$hour,jointTimes$pseudotime_shifted,method = "spearman")
fr=grpstats(jointTimes[,"pseudotime",drop=F], round(jointTimes$hour), "median")$median
fr[,1]=round(10*(fr[,1]-min(fr[,1])))
col=rainbow(max(fr)*1.2)[1:max(fr)]
col=col[fr]
# pdf("~/Downloads/realtime_vs_pseudotime.pdf",width = 4.5,height = 4.5)
vioplot::vioplot(jointTimes$pseudotime_shifted ~ round(jointTimes$hour), col=col)
legend("topleft",as.character(fr[,1]), fill=col, title="pseudotime")
# dev.off()


## Shifting real time instead of pseudotime
# jointTimes=read.csv(file=paste0(OUTPSEUDOTIME,filesep,FoV,"_matlabOut.csv"),check.names = F,stringsAsFactors = F)
# # te=cor.test(jointTimes$time_since_division,jointTimes$pseudotime_shifted,method = "spearman")
# te=cor.test(jointTimes$hour_shifted,jointTimes$pseudotime,method = "spearman")
# RES=.5
# fr=grpstats(jointTimes[,"hour",drop=F], round(jointTimes$pseudotime/RES), "median")$median
# fr[,1]=round(1+fr[,1]-min(fr[,1]))
# fr=fr[order(as.numeric(rownames(fr))),,drop=F]
# col=rainbow(max(fr)*1.3)[1:max(fr)]
# col=col[fr]
# vioplot::vioplot(jointTimes$hour_shifted ~ round(jointTimes$pseudotime/RES), col=col)
# legend("topleft",as.character(fr[,1]), fill=col, title="hour")
