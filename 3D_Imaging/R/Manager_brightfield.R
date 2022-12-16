# conda activate r_env
#setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/R")
setwd("~/Documents/Projects/PMO/MeasuringFitnessPerClone/code/3D_Imaging/R")
source("CorrectCellposeSegmentation.R")
source("assignCompartment2Nucleus.R")
source("compareCells.R")
source("generateImageMask.R")
source("Utils.R")
eps=read.table('../dbscan/eps.txt')
library(abind)
library(matlab)
library(rgl)
library(geometry)
library(magick)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rhdf5lib", force=TRUE)
BiocManager::install("rhdf5", force=TRUE)
BiocManager::install("EBImage", force=TRUE)
library("rhdf5")
library("EBImage")

##################################################################################################################
## Constants, Settings, Input and output folders:
#ROOT="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87"
ROOT="~/Documents/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87"
#ROOT="~/Documents/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/SUM-159"
setwd(ROOT)
ZSTACK_DISTANCE=0.29
EPS=5
MINPTS=4
xydim = 255
r3dDefaults$windowRect=c(0,50, 1600, 800) 
INDIR="A04_CellposeOutput"
OUTCORRECTED="A05_PostProcessCellposeOutput"
ILASTIKINPUT="H06_IlastikInput"
OUTLINKED="A06_multiSignals_Linked"
OUTSTATS="A07_LinkedSignals_Stats"
dir.create(OUTCORRECTED)
dir.create(OUTLINKED)
dir.create(OUTSTATS)
dir.create(ILASTIKINPUT)

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

assignCompartment2Nucleus<-function(other_center_coord, nuc_center_coord,OUTD){
  dirCreate(OUTD, recursive = T, permission = "a+w")
  ZSTACK_DISTANCE=0.29
  other_coord=read.csv(file=other_center_coord,check.names = F,stringsAsFactors = F)
  nuc_coord=read.csv(file=nuc_center_coord,check.names = F,stringsAsFactors = F)
  other_coord$z=other_coord$z*ZSTACK_DISTANCE
  nuc_coord$z=nuc_coord$z*ZSTACK_DISTANCE
  
  ##@TODO: should not be necessary!
  other_center_coord=strsplit(other_center_coord,"_Cells_")[[1]][1]
  nuc_center_coord=strsplit(nuc_center_coord,"_Cells_")[[1]][1]
  
  d=flexclust::dist2(other_coord, nuc_coord)
  ##Assign mitochondria to closest nucleus
  i_nuc=apply(d,1,which.min)
  ## Aggregate mitochondria all_cells_coordinates for each nucleus 
  ## and name output file according to nucleus ID
  for(i in unique(i_nuc)){
    i_mito=which(i==i_nuc)
    dm=c()
    for(j in i_mito){
      f_mito=list.files("../All_Cells_coordinates",pattern =paste0(other_center_coord,"_cell_",j,"_"),full.names = T )
      ## @TODO: should not be necessary
      if(isempty(f_mito)){
        next
      }
      coord=read.csv(file=f_mito,check.names = F,stringsAsFactors = F)
      coord[,other_center_coord]=j
      dm=rbind(dm,coord)
    }
    #### WRITE LINKED COMPARTMENTS:
    ## Copy other compartment file (e.g. nucleus)
    f_nuc=list.files("../All_Cells_coordinates",pattern =paste0(nuc_center_coord,"_cell_",i,"_"),full.names = T )
    ## @TODO: should not be necessary
    if(isempty(f_nuc) || is.null(dm)){
      next
    }
    file.copy(f_nuc,paste0(OUTD,nuc_center_coord,"_cell_",i,"_coordinates.csv"))
    ## Also write mitochondria file
    write.csv(dm,paste0(OUTD,other_center_coord,"_cell_",i,"_coordinates.csv"), row.names = F, quote = F)
  }
}

# setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87/I04_PostProcessCellposeOutput/Cells_center_coordinates")
# OUTD="../../I05_multiOrganelles_Linked/"
# nuc_center_coord="nucleus.p_Cells_Centers.csv";
# other_center_coord="mito.p_Cells_Centers.csv"
# assignCompartment2Nucleus(other_center_coord, nuc_center_coord, OUTD)
# other_center_coord="cytoplasm.p_Cells_Centers.csv"
# assignCompartment2Nucleus(other_center_coord, nuc_center_coord, OUTD)
# 
# f=list.files(OUTD,full.names = T)
# signals_per_id=plyr::count(sapply(strsplit(f,"_"), function(x) x[length(x)-1]))
# ## keep only cells with all three signals:
# toRM=signals_per_id$x[signals_per_id$freq<3]
# for(x in toRM){
#   y=list.files(OUTD,full.names = T,pattern = paste0("_",x,"_"))
#   file.remove(y)
# }

###############################################
######Allen model performance evaluation#######
###############################################

## Input and output:
FoFs=paste0("FoF",1:7,"001001_221018_brightfield")
signals=list(nucleus.p="nucleus.p_Cells_Centers.csv",mito.p="mito.p_Cells_Centers.csv")
stats=list()
mfrow3d(nr = 1, nc = 4, sharedMouse = TRUE)  
rawimges <- images <- list()
for(FoF in FoFs){
  unlink(paste0(OUTCORRECTED,filesep,FoF),recursive=T)
  
  ###############################################
  ###### Correcting Cellpose Segmentation #######
  ###############################################
  ## First correct segmentation output
  # correctSegmentations(FoF, signals, eps)
  CorrectCellposeSegmentation(FoF,signal=names(signals)[1],INDIR,OUTCORRECTED,doplot=F,eps=EPS,minPts=MINPTS,IMPORTALLORGANELLES=T)
  # rgl.snapshot("~/Downloads/Brightfield_Timeseries.png")
  images[[FoF]]=generateImageMask(FoF, INDIR=OUTCORRECTED, OUTDIR=OUTCORRECTED,root = ROOT, xydim = xydim)
  img=bioimagetools::readTIF(paste0(INDIR,filesep,FoF,filesep,names(signals),".tif"))
  img=img[fliplr(1:nrow(img)),,1:dim(images[[FoF]])[3]]
  img=EBImage::rotate(img,-90)
  rawimges[[FoF]]=resize4Ilastik(img, xydim = xydim)
  
  
  ## Next link each predicted nucleus to its closest target nucleus
  OUTLINKED_=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)
  setwd(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates"))
  # assignCompartment2Nucleus(signals$nucleus.p, signals$nucleus.t, OUTLINKED_)
  assignCompartment2Nucleus(signals$mito.p, signals$nucleus.p, OUTLINKED_)
  setwd(ROOT)
  
  ## Compare each predicted to its linked target nucleus
  #stats[[FoF]]=compareCells(signals$nucleus.t, signals$nucleus.p, OUTLINKED_)
}
#save(file="~/Documents/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87/H06_IlastikInput/221018/20055_stats.RObj","stats")

## save as h5 for Ilastik
saveAsH5<-function(images, H5OUT, binary=F, dotrim=F){
  h5=do.call(abind,c(images,along=4))
  # h5=aperm(h5,c(4,1,2,3))
  if(binary){
    h5[h5>0]=1
  }
  file.remove(H5OUT)
  h5createFile(H5OUT)  
  if(dotrim){
    h5=255*h5
    h5=h5-min(h5)
  }
  h5=round(h5)
  h5createDataset(H5OUT, dataset = fileparts(H5OUT)$name, dims=dim(h5), H5type="H5T_NATIVE_UINT32")
  h5write(h5, file = H5OUT, fileparts(H5OUT)$name)
  h5closeAll()
  return(h5)
}
tmp=gsub(substr(FoF,1,4), "FoFX",FoF)
h5=saveAsH5(rawimges,paste0(ILASTIKINPUT,filesep,tmp,".h5"),dotrim=T)
h5=saveAsH5(images,paste0(ILASTIKINPUT,filesep,tmp,"_mask.h5"))

#################################
###### Linking organelles #######
#################################
## Input and output:
FoFs=paste0("FoF",1:7,"002002_221018_brightfield")
signals=list(nucleus.p="nucleus.p_Cells_Centers.csv",mito.p="mito.p_Cells_Centers.csv")
#signals=list(nucleus.p="nucleus.p_Cells_Centers.csv",mito.p="mito.p_Cells_Centers.csv")#,cytoplasm.t="cytoplasm.t_Cells_Centers.csv")
signals_per_id=list()
for(FoF in FoFs){
  
  #OUTLINKED_=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)
  #setwd(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates"))
  # ## link each predicted nucleus to its closest target nucleus
  #assignCompartment2Nucleus(signals$mito.p, signals$nucleus.p, OUTLINKED_)
  # # assignCompartment2Nucleus(signals$cytoplasm.t, signals$nucleus.p, OUTLINKED_)
  # setwd(ROOT)
  
  ## keep only cells with all three signals:
  f=list.files(OUTLINKED_,full.names = T)
  signals_per_id_=plyr::count(sapply(strsplit(f,"_"), function(x) x[length(x)-1]))
  toRM=signals_per_id_$x[signals_per_id_$freq<length(signals)]
  for(x in toRM){
    y=list.files(OUTLINKED_,full.names = T,pattern = paste0("_",x,"_"))
    file.remove(y)
  }
  signals_per_id[[FoF]]=signals_per_id_[!signals_per_id_$x %in% toRM,]
}

###########################
## Calculate image stats ##
###########################
for(FoF in names(signals_per_id)){
  print(FoF)
  OUTLINKED_=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)
  
  ## Read in linked organelles
  coord_=readOrganelleCoordinates(signals_per_id[[FoF]], names(signals), OUTLINKED_)
  save(file=paste0("~/Downloads/stats_",FoF,".RObj"),"coord_")
  coord=coord_
  
  ## Gather stats
  cells=unique(coord_$id)
  thesignals=unique(coord_$signal)
  imgStats=as.data.frame(matrix(NA,length(cells),3+4*length(thesignals)))
  rownames(imgStats)=as.character(cells)
  colnames(imgStats)=c(sapply(c("vol_","area_","pixels_","count_"),paste0,thesignals),"x","y","z")
  for(id in cells){
    # print(paste("cell",id))
    a=coord_[coord_$id==id,]
    imgStats[as.character(id),c("x","y","z")]=apply(a[a$signal=="nucleus.p",c("x","y","z")],2,median)
    for(signal in thesignals){
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
  for(organelle in thesignals){
    imgStats[,paste0("pixels_per_volume_",organelle)]=imgStats[,paste0("pixels_",organelle)]/imgStats[,paste0("vol_",organelle)]
  }
  if("cytoplasm.t" %in% names(signals)){
    imgStats$nuc_to_mito_plus_cyto = imgStats$vol_nucleus.p/(imgStats$vol_mito.p+imgStats$vol_cytoplasm.t)
    imgStats$nuc_to_cyto = imgStats$vol_nucleus.p/imgStats$vol_cytoplasm.t
    imgStats$cyto_to_mito = imgStats$vol_cytoplasm.t/imgStats$vol_mito.p
  }
  imgStats$nuc_to_mito = imgStats$vol_nucleus.p/imgStats$vol_mito.p
  imgStats$nuc_vol_to_area=imgStats$vol_nucleus.p/imgStats$area_nucleus.p
  imgStats$z=imgStats$z/ZSTACK_DISTANCE
  write.table(imgStats,file=paste0(OUTSTATS,filesep,FoF,"_stats.txt"),sep="\t",quote=F,row.names = T)
  
  
  ## Visualize stats
  imgStats=read.table(paste0(OUTSTATS,filesep,FoF,"_stats.txt"),header = T,sep="\t")
  ## Z-score
  imgStats=apply(imgStats, 2, function(x) (x - mean(x,na.rm=T))/sd(x,na.rm=T))
  tmp=as.matrix(imgStats)
  tmp[!is.finite(tmp)]=NA
  hm = gplots::heatmap.2(tmp,trace = "none", margins = c(13, 6), symm = F)
}

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