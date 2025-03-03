# conda activate r_env
# setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/R")
setwd("/Volumes/Expansion/Collaboration/Moffitt_Noemi/BioinformaticsPaper/PMO_MeasuringFitnessperClone/3D_Imaging/R")
source("CorrectCellposeSegmentation.R")
source("assignCompartment2Nucleus.R")
source("compareCells.R")
# source("clusterMito.R")
source("generateImageMask.R")
source("Utils.R")
source("visualizeSingleCells.R")
eps=read.table('../dbscan/eps.txt')
library(abind)
library(matlab)
library(rgl)
library(ijtiff)
library(geometry)
library(habtools)
library(flexclust)
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")


## Constants, Settings, Input and output folders:
#ROOT="/Volumes/Expansion/Collaboration/Moffitt_Noemi/BioinformaticsPaper/data/NCI-N87"
ROOT="/Volumes/Expansion/Collaboration/Moffitt_Noemi/BioinformaticsPaper/data/NCI-N87"
# ROOT="~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCLs/3Dbrightfield/SUM-159"
setwd(ROOT)
EPS=10; #6
MINPTS=3; #4
xydim = 255
pixelsize_xy = 0.232 # um 
z_interval = 0.29 #  um 
xyz=c("x","y","z")
# r3dDefaults$windowRect=c(0,50, 1600, 800) 
r3dDefaults$windowRect=c(0,50, 500,500)
INDIR="A04_CellposeOutput"
OUTCORRECTED="A05_PostProcessCellposeOutput"
dirCreate(OUTCORRECTED, permission = "a+w")
# ILASTIKINPUT="G06_IlastikInput"
# dirCreate(ILASTIKINPUT, permission = "a+w")
OUTLINKED="A06_multiSignals_Linked"
OUTSTATS="A07_LinkedSignals_Stats"
dirCreate(OUTLINKED, permission = "a+w")
dirCreate(OUTSTATS, permission = "a+w")
xmlfiles=list.files('A01_rawData/',pattern=".xml",full.names=T)


## Local helper functions
correctSegmentations<-function(FoF, signals, eps){
  sapply(names(signals), function(x) CorrectCellposeSegmentation(FoF,signal=x,INDIR,OUTCORRECTED,doplot=F,eps=eps[FoF,x]))
}
readOrganelleCoordinates<-function(signals_per_id, signals, IN){
  coord=NULL;
  for(cell in signals_per_id$x){
    for(s in signals){
      x=paste0(s,"_cell_",cell,"_coordinates.csv")
      a=read.csv(paste0(IN,filesep,x))
      ## pixel to um conversion: @TODO
      a$z=a$z*z_interval
      a$x=a$x*pixelsize_xy
      a$y=a$y*pixelsize_xy
      # id=strsplit(fileparts(x)$name,"_")[[1]]
      a$organelle = a[,ncol(a)] 
      a$intensity = a$signal
      a$signal=s
      # id=id[length(id)-1]
      a$id=cell
      ## add missing columns
      if(!is.null(coord)){
        for(mc in setdiff(colnames(a), colnames(coord))){
          coord[,mc]=NA
        }
        tmp=matrix(NA,nrow(a), ncol(coord))
        colnames(tmp) = colnames(coord)
        coord=rbind(tmp, coord)
        coord[1:nrow(a),colnames(a)] = a
      }else{
        coord=a
      }
      # a[,c("y", "x", "z", "organelle", "id","signal")]
    }
  }
  coord$id=as.numeric(coord$id)
  return(coord)
}


#REGEX="240918_fluorescent.nucleus"
REGEX="231005_fluorescent.nucleus"
if (grepl("^240918", REGEX)){
  MINGREEN = 250
  MINRED = 600
  print(REGEX)
  print(sprintf("MINGREEN is %f MINRED is %f", MINGREEN,MINRED))
}else if(grepl("^231005", REGEX)){
  MINGREEN = 700
  MINRED = 500
  print(REGEX)
  print(sprintf("MINGREEN is %f MINRED is %f", MINGREEN,MINRED))
}


##########################
### read fucci results ###
FUCCIDIR=paste0("/Volumes/Expansion/Collaboration/Moffitt_Noemi/BioinformaticsPaper/data")
fucci=read.csv(list.files(FUCCIDIR,recursive = T, pattern="object_231005.csv", full.names = T)) #object_240918.csv
fucci$FileName_bright=gsub(".ome.tif","",gsub("stk_0001_","",fucci$FileName_bright))
colnames(fucci)=gsub("Location_Center_","", colnames(fucci))
ii=which(colnames(fucci) %in% toupper(xyz))
colnames(fucci)[ii]=tolower(colnames(fucci)[ii])
fucci$FileName_bright=gsub("_ch1","",fucci$FileName_bright)

###############################################
######Allen model performance evaluation#######
###############################################

## Input and output:
# FoFs=list.files(INDIR, pattern=REGEX)
FoFs=list.files(INDIR, pattern=REGEX)
signals=list(nucleus.p="nucleus.p_Cells_Centers.csv");#,mito.p="mito.p_Cells_Centers.csv", cytoplasm.p="cytoplasm.p_Cells_Centers.csv")
# signals=list(nucleus.p="nucleus.p_Cells_Centers.csv"); #nucleus.t="nucleus.t_Cells_Centers.csv",
# signals=list(nucleus.t="nucleus.t_Cells_Centers.csv"); 
stats=list()
mfrow3d(nr = 1, nc = 4, sharedMouse = TRUE)  
rawimges <- images <- list()
ncells=list()
for(FoF in FoFs){
  print(FoF)
  setwd(ROOT)
  unlink(paste0(OUTCORRECTED,filesep,FoF),recursive=T)
  
  ###############################################
  ###### Correcting Cellpose Segmentation #######
  ###############################################
  CorrectCellposeSegmentation(FoF,signal=names(signals)[1],INDIR,OUTCORRECTED,doplot=0,eps=EPS,minPts=MINPTS,IMPORTALLORGANELLES=F)
  OUTLINKED_=paste0(getwd(),filesep,OUTCORRECTED,filesep,FoF,filesep,"All_Cells_coordinates")
  ncells[[FoF]]= length(list.files(paste0(OUTCORRECTED,filesep,FoF,filesep,"All_Cells_coordinates"),pattern = "nucleus"))
  # # ## For live-cell tracking:
  # # # # CorrectCellposeSegmentation(FoF,signal=names(signals),INDIR,OUTCORRECTED,doplot=F,eps=EPS,minPts=MINPTS,IMPORTALLORGANELLES=F)
  # # # # rgl.snapshot("~/Downloads/Brightfield_Timeseries.png")
  # # # images[[FoF]]=generateImageMask(FoF, INDIR=OUTCORRECTED, OUTDIR=OUTCORRECTED,root = ROOT, xydim = xydim)
  # images[[FoF]]=try(generateImageMask(FoF, INDIR=OUTCORRECTED, OUTDIR=OUTCORRECTED,root = ROOT, signal = "nucleus.p"))
  # # # img=bioimagetools::readTIF(paste0(INDIR,filesep,FoF,filesep,names(signals),".tif"))
  # # # img=img[fliplr(1:nrow(img)),,1:dim(images[[FoF]])[3]]
  # # # img=EBImage::rotate(img,-90)
  # # # rawimges[[FoF]]=resize4Ilastik(img, xydim = xydim)
  # 
  
  # ## Cluster Organelles
  # setwd(ROOT)
  # OUTSEG_=paste0(getwd(),filesep,OUTCORRECTED,filesep,FoF,filesep)
  # ##@TODO rm all mito and cyto files
  MITOTIF=paste0(ROOT,filesep,INDIR,filesep,FoF,"/mito.p.tif")
  # clusterMito(MITOTIF, OUTSEG_)
  CYTOTIF=paste0(ROOT,filesep,INDIR,filesep,FoF,"/cytoplasm.p.tif")
  # clusterMito(CYTOTIF, OUTSEG_)
  
  # ## For linking organelles or to link each predicted nucleus to its closest target nucleus
  setwd(ROOT)
  OUTLINKED_=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)
  unlink(OUTLINKED_,recursive=T)
  setwd(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates"))
  # # # assignCompartment2Nucleus(signals$nucleus.p, signals$nucleus.t, OUTLINKED_)
  # assignCompartment2Nucleus(signals$mito.p, signals$nucleus.p, OUTLINKED_, save_cell_gif=F)
  # assignCompartment2Nucleus(signals$cytoplasm.p, signals$nucleus.p, OUTLINKED_, save_cell_gif=F)
  assignCompartment2Nucleus(MITOTIF, CYTOTIF,OUTLINKED_, signals$nucleus.p, save_cell_gif=F)
  setwd(ROOT)
  
  
  ## Save FUCCI intensities
  ## Use Nucleus coordinates (from cellpose postprocess) to access intensities in ch00 and ch02 of FUCCI image:
  fucci_coord=list()
  for(ch in c("ch00","ch02")){
    TIF=list.files(paste0("A01_rawData",filesep,FoF), pattern=ch, full.names = T)
    ch_=sapply(TIF, function(x) bioimagetools::readTIF(x), simplify = F)
    fucci_coord[[ch]]=do.call(abind, ch_)
  }
  f=list.files(OUTLINKED_, pattern="nucleus.p_cell_",full.names = T)
  ## Save cellprofiler output with matching nucleus ID:
  fucci_=fucci[fucci$FileName_bright==FoF,]
  fucci_$ID=NA
  colnames(fucci_) = gsub("fluor_1","green",colnames(fucci_)) ## fluor_1 = green
  colnames(fucci_) = gsub("fluor_2","red",colnames(fucci_)) ## fluor_2 = red
  ## classify cell cycle phase
  fucci_$cellCycle = 2
  fucci_$cellCycle[fucci_$Intensity_IntegratedIntensity_green>MINGREEN & fucci_$Intensity_IntegratedIntensity_red<MINRED] = 1
  fucci_$cellCycle[fucci_$Intensity_IntegratedIntensity_green<MINGREEN & fucci_$Intensity_IntegratedIntensity_red>MINRED] = 3
  fucci_$cellCycle[fucci_$Intensity_IntegratedIntensity_green>MINGREEN & fucci_$Intensity_IntegratedIntensity_red>MINRED] = 4
  for(x in f){
    coord=read.csv(file=x,check.names = F,stringsAsFactors = F)
    coord_ = coord
    coord_[,xyz] = coord_[,xyz]+1
    coord$fucci_ch00=apply(coord_,1, function(p) fucci_coord$ch00[p["y"],p["x"],p["z"]])
    coord$fucci_ch02=apply(coord_,1, function(p) fucci_coord$ch02[as.numeric(p["y"]),as.numeric(p["x"]),as.numeric(p["z"])] )
     ## Assign FUCCI ID:
    d=dist2(t(as.matrix(apply(coord_[,xyz],2,mean))), fucci_[,xyz])
    fucci_$ID[which.min(d)] = as.numeric(strsplit(fileparts(x)$name,"_")[[1]][3])
    
    coord$cellCycle = fucci_$cellCycle[which.min(d)]
    ## Save output -- overwrite
    write.csv(coord, file=x,quote = F,row.names = F)
  }
  fucci_=fucci_[!is.na(fucci_$ID),]
  rownames(fucci_)=fucci_$ID
  write.table(fucci_, file=paste0(OUTLINKED_,filesep,FoF,"_fucci.txt"), quote = F,row.names = F, sep="\t")
  
  # ## Visualize cells
  # cells=unique(sapply(strsplit(list.files(OUTLINKED_,pattern = "nucleus.p"),"_"),"[[",3))
  # tmp=sapply(cells, function(i) visualizeSingleCells(i, signals$mito.p, signals$nucleus.p, OUTLINKED_))
  
  # ## Compare each predicted to its linked target nucleus
  # stats[[FoF]]=compareCells(signals$nucleus.t, signals$nucleus.p, OUTLINKED_)
}


## keep only cells with all three signals:
signals=list(nucleus.p="nucleus.p_Cells_Centers.csv",mito.p="mito.p_Cells_Centers.csv",cytoplasm.p="cytoplasm.p_Cells_Centers.csv")
signals_per_id=list()
## Input and output:
for(FoF in FoFs){
  OUTLINKED_=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)
  
  f=list.files(OUTLINKED_,full.names = T, pattern = ".csv")
  signals_per_id_=plyr::count(sapply(strsplit(f,"_"), function(x) x[length(x)-1]))
  toRM=signals_per_id_$x[signals_per_id_$freq<length(signals)]
  for(x in toRM){
    y=list.files(OUTLINKED_,full.names = T,pattern = paste0("_",x,"_"))
    file.remove(y)
  }
  signals_per_id[[FoF]]=signals_per_id_[!signals_per_id_$x %in% toRM,]
}
