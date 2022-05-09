# conda activate r_env
setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/R")
source("CorrectCellposeSegmentation.R")
eps=read.table('../dbscan/eps1.txt')
library(matlab)
library(rgl)
library(geometry)

## Constants, Settings, Input and output folders:
ROOT="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87"
setwd(ROOT)
ZSTACK_DISTANCE=0.29
r3dDefaults$windowRect=c(0,50, 800, 800) 
FoFs="FoF12_211110_fluorescent.nucleus"
EPS=eps[1,1]
INDIR="A04_CellposeOutput"
#OUTCORRECTED=paste0("C05_TestDBSCANsettings",filesep,FoFs,filesep,EPS)
OUTCORRECTED=paste0("C05_TestDBSCANsettings",filesep,FoFs,filesep,EPS)
OUTSTATS=paste0("C06_Stats",filesep,FoFs)
dir.create(OUTCORRECTED,recursive=T)
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
###### Correcting Cellpose Segmentation #######
###############################################
## Input and output:
signals=list(nucleus.t="nucleus.t_Cells_Centers.csv")
stats=list()
for(FoF in FoFs){
  ## First correct segmentation output
  correctSegmentations(FoF, signals, eps)
}

OUTCORRECTED=paste0("C05_TestDBSCANsettings",filesep,FoFs,filesep,EPS)
#filesep,OUTCORRECTED,filesep,FoF,filesep, 
OUTCORRECTED_=paste0(getwd(),filesep,OUTCORRECTED,filesep,"All_Cells_coordinates")

## Keep only cells with one signal:
f=list.files(OUTCORRECTED_,full.names = T, pattern="nucleus.t")
signals_per_id=plyr::count(sapply(strsplit(f,"_"), function(x) x[length(x)-1]))

## Read in one organelle
coord_=readOrganelleCoordinates(signals_per_id, names(signals), OUTCORRECTED_)
coord=coord_

## Gather stats
cells=unique(coord_$id)
signals=unique(coord_$signal)
imgStats=as.data.frame(matrix(NA,length(cells),4*length(signals)))
rownames(imgStats)=as.character(cells)
colnames(imgStats)=c(sapply(c("vol_","area_","pixels_","count_"),paste0,signals))
for(id in cells){
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

# Add unique identifier columns for eps
imgStats$eps = EPS

# Save stats in csv
write.csv(imgStats,file=paste0(OUTSTATS,filesep,EPS,"eps_stats.csv"), row.names = FALSE)