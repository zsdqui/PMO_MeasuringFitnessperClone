compareCells<-function(signal1, signal2,OUTD){
  library(matlab)
  COLS=c("x","y","z")
  ZOOMOUTFACTOR=1
  types=sapply(strsplit(c(signal1,signal2),"_Cells"),"[[",1)
  
  ids=list.files(OUTD,pattern=types[1])
  ids=sapply(strsplit(ids,"_"),"[[",3)
  
  stats=as.data.frame(matrix(NA,length(ids),4))
  rownames(stats)=ids
  colnames(stats)=c(paste0(types,"_NumPixels"),paste0(types,"_IntersectingPixels"))
  
  for(id in ids){
    pair=list.files(OUTD,pattern=paste0("_",id,"_"),full.names = T)  
    coord=lapply(pair, read.csv)
    names(coord)=sapply(pair,function(x) strsplit(fileparts(x)$name,"_")[[1]][1])
    
    i=dplyr::inner_join(round(coord[[1]][,COLS]/ZOOMOUTFACTOR),round(coord[[2]][,COLS]/ZOOMOUTFACTOR))
    stats[id,paste0(names(coord),"_IntersectingPixels")]=nrow(i)/sapply(coord,nrow)
    stats[id,paste0(names(coord),"_NumPixels")]=sapply(coord,nrow)
  }
  return(stats)
}