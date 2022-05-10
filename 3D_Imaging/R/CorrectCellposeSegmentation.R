CorrectCellposeSegmentation<- function(ID,signal,INDIR,OUTDIR,doplot=F, eps=2.5, minPts = 2, IMPORTALLORGANELLES=T, MINZSLICES=0.4*70){
  library(matlab)
  library(geometry)
  library(misc3d)
  library(rgl)
  r3dDefaults$windowRect=c(0,50, 800, 800) 
  
  ## Create output folders
  sapply(c("Cells_center_coordinates","All_Cells_coordinates"), function(x) dir.create(paste0(OUTDIR,filesep,ID,filesep,x), recursive = T  ))
  
  ############################################################
  ## Correct segmentation: merge ids belonging to same cell ##
  coord_=read.csv(paste0("./",INDIR,filesep,ID, "/Cells_center_coordinates/", signal, "_Cells_Centers.csv"))
  o=dbscan::dbscan(coord_[,c("x","y")],eps = eps, minPts = minPts)
  Sys.sleep(10)
  file.remove("Rplots.pdf")
  coord_$id=o$cluster
  fr=plyr::count(o$cluster)
  colnames(fr)[1]="newCellID"
  fr=fr[order(fr$freq,decreasing = T),]
  fr=fr[fr$newCellID!=0,]; ##remove noise
  
  fr$datapoints=NA
  newCellCoord=list()
  allNewCellIDs = fr$newCellID
  for(newCell in fr$newCellID){
    ##Gather all coordinates associated with new cell
    oldCells=which(newCell==coord_$id)
    tmp=c()
    for(oldid in oldCells){
      dm_=read.csv(paste0("./",INDIR,filesep,ID, "/All_Cells_coordinates/", signal, "_cell_",oldid,"_coordinates.csv"))
      dm_$oldid=oldid
      tmp=rbind(tmp,dm_)
    }
    ## Filter by zstack representation
    if(length(round(unique(tmp$z),1)) < MINZSLICES){
      allNewCellIDs =setdiff(allNewCellIDs, newCell)
      next
    }
    newCellCoord[[as.character(newCell)]]=tmp
    ## Record # of coordinates:
    fr$datapoints[fr$newCellID==newCell]=nrow(tmp)
  }
  if(isempty(allNewCellIDs)){
    print("None of the cells have sufficient zstack representation. No cells saved")
    return()
  }
  fr = fr[fr$newCellID %in% allNewCellIDs, ]
  ## Stats
  centroids=t(sapply(newCellCoord, function(x) apply(x[,c("x","y","z")],2,mean)))
  colnames(centroids)=c("x","y","z")
  fr=fr[order(fr$datapoints),]
  hist(fr$freq)
  ## Plot coordinates for cells of interest
  if(doplot){
    o2=dbscan::dbscan(centroids[,c("x","y")],eps = eps*8, minPts = 2)
    coi=rownames(centroids)[o2$cluster==3]
    plot(centroids[,"x"],-centroids[,"y"],pch=20,col=1+(rownames(centroids) %in% coi))
    # coi=fr$newCellID[seq(1,nrow(fr),by=20)]
    rgl::close3d()
    for(cell in coi){
      a=newCellCoord[[as.character(cell)]]
      hull=Plot_ConcaveHull(a$x, -a$y, a$z, lcolor =which(cell==coi), alpha=0.5,add = T)
    }
    view3d( theta = 0, phi = 0, fov=0)
  }
  ## Save output
  tmp=sapply(1:length(newCellCoord), function(id) write.csv(newCellCoord[[id]],file=paste0("./",OUTDIR,filesep,ID, "/All_Cells_coordinates/", signal, "_cell_",id,"_coordinates.csv"), row.names = F, quote = F))
  write.csv(centroids,file=paste0("./",OUTDIR,filesep,ID, "/Cells_center_coordinates/", signal, "_Cells_Centers.csv"),quote = F, row.names = F)
  print(paste(nrow(coord_),"cell IDs merged into",length(newCellCoord),"unique IDs"))
  
  ## Also copy mitochondria and cytoplasm (if they exist):
  if(IMPORTALLORGANELLES){
    for(other in c("mito","cyto")){
      for(res in c("Cells_center_coordinates","All_Cells_coordinates")){
        sapply(list.files(paste0(INDIR,filesep,ID,"/",res,"/"),pattern = other,full.names = T), function(x) file.copy(x,paste0(OUTDIR,filesep,ID, "/",res,"/",fileparts(x)$name,".csv")))
      }
    }
  }
}
