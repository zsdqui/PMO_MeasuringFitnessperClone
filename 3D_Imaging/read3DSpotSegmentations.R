setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/3Dbrightfield/allencell")
library(matlab)
library(geometry)
ZSTACK_DISTANCE=0.9
f=list.files("G03_CellposeOutput",full.names = T)

## Gather stats
imgStats=as.data.frame(matrix(NA,length(f),3))
rownames(imgStats)=sapply(f,function(x) fileparts(x)$name)
colnames(imgStats)=c("volume","area","id")
coord=c();
for(x in f){
  a=read.csv(x)
  a$z=a$z*ZSTACK_DISTANCE
  hull <- convhulln(a, options = "FA")
  imgStats[fileparts(x)$name,c("id","volume","area")]=c(id,hull$vol,hull$area)
  
  id=strsplit(fileparts(x)$name,"_")[[1]]
  id=id[length(id)-1]
  a$id=id
  coord=rbind(coord,a)
}
write.table(imgStats,file="G04_segmentationStats/0_prediction_c0.model_stats.txt",sep="\t")
# apply(coord[,1:3],2,quantile,c(0,1))

###########
## Plots ##
## plot all cells
realCells=imgStats$id[imgStats$volume>=100]
coord_=coord[coord$id %in% realCells,]
o=kmeans(coord_[,c("x","y","z")],centers = 12)
tmp = quantile(coord_$z,c(0,1))
space=tmp[2]-tmp[1]
rgl.close()
gl::plot3d(coord_$x, coord_$y, coord_$z, pch=20, zlim=c(tmp[1]-space/2, tmp[2]+space/2), size=2, axes=F, xlab="",ylab="", zlab="",col=o$cluster+1)
plot(1);legend("topleft",as.character(unique(o$cluster)),fill = unique(o$cluster)+1)
rgl.close()

## plot concave hull for one cell cluster
fr=plyr::count(o$cluster)
fr=fr[order(fr$freq),]
rgl::close3d()
for(id in unique(coord_$id[o$cluster==6])){
  a=coord_[coord_$id==id,]
  hull=Plot_ConcaveHull(a[,1], a[,2], a[,3], lcolor ="cyan", alpha=0.95,add = T)
}
rgl::movie3d(
  movie="CellPose3D_output",
  rgl::spin3d( axis = c(1, 1, 1), rpm = 3),
  duration = 10,
  dir = "~/Downloads/",
  type = "gif",
  clean = TRUE
)