# conda activate r_env
setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/R")
source("CorrectCellposeSegmentation.R")
source("generateImageMask.R")
eps=read.table('../dbscan/eps1.txt')
library(matlab)
library(rgl)
library(geometry)


## Constants, Settings, Input and output folders:
ROOT="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87"
setwd(ROOT)
# ZSTACK_DISTANCE=0.29
PIXEL2UM = list(X=246.03/1024, Y=246.03/1024, Z=21/70)
MINZSLICES = 0.4*70
r3dDefaults$windowRect=c(0,50, 800, 800) 
FoF="FoF12_211110_fluorescent.nucleus"
INDIR="A04_CellposeOutput"

## Local helper functions
correctSegmentations<-function(FoF, signals, eps, minPts){
  sapply(names(signals), function(x) CorrectCellposeSegmentation(FoF,signal=x,INDIR,OUTCORRECTED,doplot=F,eps=eps,minPts=minPts,IMPORTALLORGANELLES=F))
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
  # coord$z=coord$z*ZSTACK_DISTANCE
  return(coord)
}
visualize3Dsegmentations<-function(coord_,prefix){
  ## Visualize 3D segmentations
  for(theta in seq(0,90,by=30)){ #seq(30,90,by=30)
    for(phi in seq(0,90,by=30)){
      tiff(paste0(prefix,"_3D_",theta,"_",phi,".tif"),width = 1080, height = 1080)
      plot3D::scatter3D(coord_$x, coord_$y, coord_$z, size=2, axes=F, xlab="",ylab="", zlab="",colvar=coord_$id,transparency=0.2,theta=theta,phi=phi)
      dev.off()
    }
  }
  cmd=paste0("convert -delay 20 -loop 0 ",prefix,"_3D_*.tif ",prefix,"_3D.gif ")
  system(cmd)
  
  ## Clean up:
  sapply(list.files(fileparts(prefix)$pathstr, pattern=".tif", full.names=T), file.remove)
  # rgl::plot3d(coord_$x, coord_$y, coord_$z, pch3d=20, size=2, axes=F, zlim=zlim, xlab="",ylab="", zlab="",col=coord_$id+1,alpha=0.04)
  # rgl::rgl.snapshot(paste0(OUTSTATS,filesep,FoF,"_3D.png"), fmt = 'png')
  # ## Save as gif
  # rgl::movie3d(
  #   movie=paste0(FoF,"_3D"),
  #   rgl::spin3d( axis = c(1, 1, 1), rpm = 8),
  #   duration = 1,
  #   dir = OUTSTATS,
  #   type = "gif",
  #   clean = TRUE
  # )
}



eps[,2:10]=seq(1,18,by=2)
minpts=t(as.data.frame(seq(2,4, by=1)))
rownames(minpts) = rownames(eps)
for(EPS in sort(eps[FoF,],decreasing = F)){
  for(MINPTS in minpts[FoF,]){
    OUTCORRECTED=paste0("D05_TestDBSCANsettings",filesep,"EPS",EPS,"_MINPTS",MINPTS)
    OUTSTATS=paste0("D06_Stats",filesep,"EPS",EPS,"_MINPTS",MINPTS)
    dirCreate(OUTCORRECTED,recursive=T,permission = "a+w")
    dirCreate(OUTSTATS,recursive=T,permission = "a+w")
    unlink(paste0(OUTCORRECTED,filesep,FoF),recursive=T)
    unlink(paste0(OUTSTATS,filesep,FoF),recursive=T)
    
    ###############################################
    ###### Correcting Cellpose Segmentation #######
    ###############################################
    signals=list(nucleus.t="nucleus.t_Cells_Centers.csv")
    ## First correct segmentation output
    correctSegmentations(FoF, signals, eps=EPS, minPts=MINPTS)
    generateImageMask(FoF, INDIR=OUTCORRECTED, OUTDIR=OUTSTATS)
    
    
    ## Keep only cells with one signal:
    OUTCORRECTED_=paste0(getwd(),filesep,OUTCORRECTED,filesep,FoF,filesep, "All_Cells_coordinates", filesep)
    f=list.files(OUTCORRECTED_,full.names = T, pattern="nucleus.t")
    signals_per_id=plyr::count(sapply(strsplit(f,"_"), function(x) x[length(x)-1]))
    print(paste("Calculating stats for",OUTCORRECTED_))
    
    ## Read in one organelle
    coord_=readOrganelleCoordinates(signals_per_id, names(signals), OUTCORRECTED_)
    coord=coord_
    visualize3Dsegmentations(coord_,prefix=paste0(OUTSTATS,filesep,FoF))
    
    ## pixel2um conversion
    coord_$x = coord_$x*PIXEL2UM$X
    coord_$y = coord_$y*PIXEL2UM$Y
    coord_$z = coord_$z*PIXEL2UM$Z
    
    ## Gather stats
    cells=unique(coord_$id)
    signals=unique(coord_$signal)
    imgStats=as.data.frame(matrix(NA,length(cells),5*length(signals)))
    rownames(imgStats)=as.character(cells)
    colnames(imgStats)=c(sapply(c("vol_","area_","pixels_","zslices_","count_"),paste0,signals))
    for(id in cells){
      a=coord_[coord_$id==id,]
      
      for(signal in signals){
        stats_=list(area_=NA,vol_=NA,pixels_=NA, zslices_=NA)
        for(organelle in unique(a$organelle)){
          a_=a[a$organelle==organelle & a$signal==signal,]
          hull <- try(convhulln(a_[,c("x","y","z")], options = "FA"),silent = T)
          if(class(hull)!="try-error"){
            stats_$area_=c(stats_$area_,hull$area)
            stats_$vol_=c(stats_$vol_,hull$vol)
            stats_$pixels_=c(stats_$pixels_,nrow(a_))
          }
          stats_$zslices_=c(stats_$zslices_,length(unique(a_$z)))
        }
        imgStats[as.character(id),paste0(c(names(stats_),"count_"),signal)]=c(sapply(stats_,sum,na.rm=T),length(stats_$area_)-1)
      }
    }
    
    # Add unique identifier columns for eps
    imgStats$eps = EPS
    imgStats$minPts = MINPTS
    
    # Save stats in csv
    write.csv(imgStats,file=paste0(OUTSTATS,filesep,FoF,"_stats.csv"), row.names = FALSE)
    
    ## Generate masks for cells with sufficient z-stack representation only
    cells = rownames(imgStats)[imgStats$zslices_nucleus.t>=MINZSLICES]
    generateImageMask(FoF, INDIR=OUTCORRECTED, OUTDIR=OUTSTATS, targetcellids=cells)
  }
}