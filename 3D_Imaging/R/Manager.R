# conda activate r_env
# setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/R")
setwd("~/Projects/PMO/MeasuringFitnessPerClone/code/3D_Imaging/R")
source("CorrectCellposeSegmentation.R")
source("assignCompartment2Nucleus.R")
source("compareCells.R")
eps=read.table('../dbscan/eps.txt')
library(matlab)
library(rgl)
library(geometry)


## Constants, Settings, Input and output folders:
# ROOT="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87"
ROOT="~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87"
setwd(ROOT)
ZSTACK_DISTANCE=0.29
r3dDefaults$windowRect=c(0,50, 800, 800) 
INDIR="A04_CellposeOutput"
OUTCORRECTED="A05_PostProcessCellposeOutput"
OUTLINKED="A06_multiSignals_Linked"
OUTSTATS="A07_LinkedSignals_Stats"
dir.create(OUTCORRECTED)
dir.create(OUTLINKED)
dir.create(OUTSTATS)
## Local helper functions
correctSegmentations<-function(FoF, signals, eps){
  sapply(names(signals), function(x) CorrectCellposeSegmentation(FoF,signal=x,INDIR,OUTCORRECTED,doplot=F,eps=eps[FoF,x]))
}
readOrganelleCoordinates<-function(signals_per_id, signals, IN){
  coord=c();
  for(cell in signals_per_id$x){
    for(s in signals){
      x=paste0(s,"_cell_",cell,"_coordinates.csv")
      a=read.csv(paste0(IN,filesep,x))
      id=strsplit(fileparts(x)$name,"_")[[1]]
      a$organelle = a[,ncol(a)] 
      a$signal=s
      id=id[length(id)-1]
      a$id=id
      coord=rbind(coord,a[,c("y", "x", "z", "organelle", "id","signal")])
    }
  }
  coord$id=as.numeric(coord$id)
  coord$z=coord$z*ZSTACK_DISTANCE
  return(coord)
}


###############################################
######Allen model performance evaluation#######
###############################################

## Input and output:
FoFs=c("FoF0_211007_fluorescent.nucleus", "FoF10_210803_fluorescent.nucleus")
signals=list(nucleus.t="nucleus.t_Cells_Centers.csv", nucleus.p="nucleus.p_Cells_Centers.csv")
stats=list()
for(FoF in FoFs){
  ## First correct segmentation output
  correctSegmentations(FoF, signals, eps)
  
  ## Next link each predicted nucleus to its closest target nucleus
  OUTLINKED_=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)
  setwd(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates"))
  assignCompartment2Nucleus(signals$nucleus.p, signals$nucleus.t, OUTLINKED_)
  setwd(ROOT)
  
  ## Compare each predicted to its linked target nucleus
  stats[[FoF]]=compareCells(signals$nucleus.t, signals$nucleus.p, OUTLINKED_)
}
save(file="~/Downloads/stats.RObj","stats")

##Plot stats for first FoF
stats_=stats[[1]]
minmax=quantile(unlist(stats_[,1:2]), c(0,1))
par(mfrow=c(2,2))
plot(stats_$nucleus.t_NumPixels,stats_$nucleus.p_NumPixels,pch=20,log="xy",xlim=minmax,ylim=minmax)
hist(stats_$nucleus.t_IntersectingPixels,col="cyan")
hist(stats_$nucleus.p_IntersectingPixels,col="cyan")



#################################
###### Linking organelles #######
#################################

## Input and output:
# FoF="FoF13_220228_fluorescent.cytoplasm"
FoF="FoF16_210818_fluorescent.cytoplasm"
signals=list(nucleus.p="nucleus.p_Cells_Centers.csv",mito.p="mito.p_Cells_Centers.csv",cytoplasm.t="cytoplasm.t_Cells_Centers.csv")
OUTLINKED_=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)

## First correct segmentation output
eps=matrix(5,1,length(signals))
rownames(eps)=FoF
colnames(eps)=names(signals)
correctSegmentations(FoF, signals["nucleus.p"], eps)

## Next link each predicted nucleus to its closest target nucleus
setwd(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates"))
assignCompartment2Nucleus(signals$mito.p, signals$nucleus.p, OUTLINKED_)
assignCompartment2Nucleus(signals$cytoplasm.t, signals$nucleus.p, OUTLINKED_)
setwd(ROOT)

## keep only cells with all three signals:
f=list.files(OUTLINKED_,full.names = T)
signals_per_id=plyr::count(sapply(strsplit(f,"_"), function(x) x[length(x)-1]))
toRM=signals_per_id$x[signals_per_id$freq<length(signals)]
for(x in toRM){
  y=list.files(OUTLINKED_,full.names = T,pattern = paste0("_",x,"_"))
  file.remove(y)
}
signals_per_id=signals_per_id[!signals_per_id$x %in% toRM,]

## Read in linked organelles
coord_=readOrganelleCoordinates(signals_per_id, names(signals), OUTLINKED_)
save(file=paste0("~/Downloads/stats_",FoF,".RObj"),"coord_")
coord=coord_


##########################
## Visualize organelles ##
##########################
tmp = quantile(coord_$z,c(0,1))
space=tmp[2]-tmp[1]
zlim=c(tmp[1]-space/2, tmp[2]+space/2)
col=rainbow(length(unique(coord_$signal)))
names(col)=as.character(unique(coord_$signal))
## Color by organelle
rgl::close3d()
# rgl::plot3d(coord_$x, coord_$y, coord_$z, pch3d=20, zlim=zlim, size=2, axes=F, xlab="",ylab="", zlab="",col="white",alpha=0.4)
alpha=list(cytoplasm.p=0.01,cytoplasm.t=0.01,nucleus.t=1,nucleus.p=1,mito.t=0.1,mito.p=0.1)
for(s in names(signals)){
  X=coord_[coord_$signal==s,]
  if(s=="nucleus.p"){
    rgl::plot3d(X$x, X$y, X$z, pch3d=20, zlim=zlim, size=2, axes=F, xlab="",ylab="", zlab="",col=col[X$signal],alpha=alpha[[s]], add=T)
  }else{
    rgl::points3d(X$x, X$y, X$z, pch3d=20, col=col[X$signal],alpha=alpha[[s]], add=T)
  }
}
## Color by cell
rgl::close3d()
col=rainbow(length(unique(coord_$id)))
names(col)=as.character(unique(coord_$id))
rgl::plot3d(coord_$x, coord_$y, coord_$z, pch=20, zlim=zlim, size=2, axes=F, xlab="",ylab="", zlab="",col=col[as.character(coord_$id)],add=T)
## Save as gif
# rgl::movie3d(
#   movie=paste0("CellPose3D_output_",FoF),
#   rgl::spin3d( axis = c(1, 1, 1), rpm = 8),
#   duration = 1,
#   dir = "~/Downloads/",
#   type = "gif",
#   clean = TRUE
# )

## Gather stats
cells=unique(coord_$id)
signals=unique(coord_$signal)
imgStats=as.data.frame(matrix(NA,length(cells),4*length(signals)))
rownames(imgStats)=as.character(cells)
colnames(imgStats)=c(sapply(c("vol_","area_","pixels_","count_"),paste0,signals))
for(id in cells){
  # print(paste("cell",id))
  a=coord_[coord_$id==id,]
  for(signal in signals){
    stats_=list(area_=NA,vol_=NA,pixels_=NA)
    for(organelle in unique(a$organelle)){
      a_=a[a$organelle==organelle & a$signal==signal,]
      hull <- try(convhulln(a_[,c("x","y","z")], options = "FA"),silent = T)
      if(class(hull)!="try-error"){
        stats_$area_=c(stats_$area_,hull$area)
        stats_$vol_=c(stats_$vol_,hull$vol)
        stats_$pixels_=c(stats_$pixels_,nrow(a_))
      }
    }
    imgStats[as.character(id),paste0(c(names(stats_),"count_"),signal)]=c(sapply(stats_,sum,na.rm=T),length(stats_$area_)-1)
  }
}
## Add more stats and save 
imgStats$pixel_per_mito_avg=imgStats$pixels_mito.p/imgStats$count_mito.p
imgStats$pixel_per_volume_mito=imgStats$pixels_mito.p/imgStats$vol_mito.p
for(organelle in signals){
  imgStats[,paste0("pixels_per_volume_",organelle)]=imgStats[,paste0("pixels_",organelle)]/imgStats[,paste0("vol_",organelle)]
}
imgStats$nuc_to_mito_plus_cyto = imgStats$vol_nucleus.p/(imgStats$vol_mito.p+imgStats$vol_cytoplasm.t)
imgStats$nuc_to_cyto = imgStats$vol_nucleus.p/imgStats$vol_cytoplasm.t
imgStats$nuc_to_mito = imgStats$vol_nucleus.p/imgStats$vol_mito.p
imgStats$cyto_to_mito = imgStats$vol_cytoplasm.t/imgStats$vol_mito.p
imgStats$nuc_vol_to_area=imgStats$vol_nucleus.p/imgStats$area_nucleus.p
write.table(imgStats,file=paste0(OUTSTATS,filesep,FoF,"_stats.txt"),sep="\t",quote=F,row.names = T)


## Visualize stats
imgStats=read.table(paste0(OUTSTATS,filesep,FoF,"_stats.txt"),header = T,sep="\t")
## Z-score
imgStats=apply(imgStats, 2, function(x) (x - mean(x,na.rm=T))/sd(x,na.rm=T))
tmp=as.matrix(imgStats)
tmp[!is.finite(tmp)]=NA
hm = gplots::heatmap.2(tmp,trace = "none", margins = c(13, 6), symm = F)

