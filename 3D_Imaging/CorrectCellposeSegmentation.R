setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87/")
library(matlab)
library(geometry)
library(misc3d)
library(rgl)
r3dDefaults$windowRect=c(0,50, 800, 800) 
ZSTACK_DISTANCE=1;#0.29

############################################################
## Correct segmentation: merge ids belonging to same cell ##
coord_=read.csv("./H03_CellposeOutput/Cells_center_coordinates/nucleus.p_Cells_Centers.csv")
o=dbscan::dbscan(coord_[,c("x","y","z")],eps = 2.5,minPts = 3)
coord_$id=o$cluster
fr=plyr::count(o$cluster)
colnames(fr)[1]="newCellID"
fr=fr[order(fr$freq,decreasing = T),]
fr=fr[fr$newCellID!=0,]; ##remove noise

fr$datapoints=NA
newCellCoord=list()
for(newCell in fr$newCellID){
  ##Gather all coordinates associated with new cell
  oldCells=which(newCell==coord_$id)
  tmp=c()
  for(oldid in oldCells){
    dm_=read.csv(paste0("./H03_CellposeOutput/All_Cells_coordinates/nucleus.p_cell_",oldid,"_coordinates.csv"))
    dm_$oldid=oldid
    tmp=rbind(tmp,dm_)
  }
  newCellCoord[[as.character(newCell)]]=tmp
  ## Record # of coordinates:
  fr$datapoints[fr$newCellID==newCell]=nrow(tmp)
}
fr=fr[order(fr$datapoints),]
hist(fr$freq)
rgl::close3d()
## Plot coordinates for cells of interest
coi=fr$newCellID[seq(1,nrow(fr),by=10)]
for(cell in coi){
  a=newCellCoord[[as.character(cell)]]
  hull=Plot_ConcaveHull(a$x, a$y, a$z, lcolor =which(cell==coi), alpha=0.5,add = T)
}
## Save output
sappply(names(newCellCoord), function(id) write.csv(newCellCoord[[id]],file=paste0("./I04_PostProcessCellposeOutput/All_Cells_coordinates/nucleus.p_cell_",id,"_coordinates.csv")))

