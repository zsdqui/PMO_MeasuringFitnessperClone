library(anticlust)
setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87")
FROMILASTIK="G07_IlastikOutput"
A01="A01_rawData"
FoF="FoFX007_220523_brightfield"
XYZ=c("Center_of_the_object_0","Center_of_the_object_1","Center_of_the_object_2")

TemporalQuant_sinceANDuntil_division<-function(X){ 
  X$first_tpc <- X$last_tpc <- X$time_of_division <- NA
  ## When does a cell first appear and when does it dissapear
  tpc=grpstats(as.matrix(X$POSIXct),X$trackId,c("min","max"))
  dpc=grpstats(as.matrix(X$POSIXct),X$parentTrackId,"min")$min
  for(cell in unique(X$trackId)){
    ii=which(X$trackId==cell)
    X$first_tpc[ii]=tpc$min[as.character(cell),]
    X$last_tpc[ii]=tpc$max[as.character(cell),]
    if(as.character(cell) %in% rownames(dpc)){
      X$time_of_division[ii] = dpc[as.character(cell),]
    }
  }
  hist(X$last_tpc-X$first_tpc,main="lifetime of cell")
  ## Plots
  par(mfrow=c(2,2))
  plot(tpc$min,ylab="first_tpc")
  plot(tpc$max,ylab="last_tpc")
  plot(dpc,ylab="time_of_division")
  # ##pick two later cells and set their parent to be one of the earlier cells
  # lateCells=as.numeric(rownames(X$first_tpc)[X$first_tpc>0])
  # earlyCells=as.numeric(rownames(X$first_tpc)[X$first_tpc==0])
  # daughterCells=lateCells[1:2]
  # parentalCell=earlyCells[1]
  # X$parentTrackId[X$trackId %in% daughterCells] = parentalCell
  return(X)
}


# Ilastik settings to obtain this result:
# Appearance cost: 10
# Disappearance Cost: 10
# Transition Weight: 20
# Division Weight: 0
# Other: default


##Center_of_the_object_0 and Center_of_the_object_1 are the x and y axis respectively
##Frame = timepoint
# labelimageId: id of the object in the current time frame. Will range from 1 .. number_of_objects in the current time frame.
# lineageId: id of the linage. All cells that are children of a certain object should have the same lineageId. A value of -1 indicates a false detection.
# trackId: id of the part of the lineage tree. A value of -1 indicates a false detection. Events like appearance, disappearance and division will start a new track.
# parentTrackId: id of the track id of the parent object track (in the previous time frame). If an object divides, there will be two new tracks belonging to the same lineage.
# X=read.csv("~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/LiveCellImaging/A01_stitchedwells-stiched_well_CSV-Table.csv")
X=read.csv(paste0(FROMILASTIK,filesep,FoF,"_CSV-Table.csv"))
print(paste("Found",nrow(X[X$trackId>=0,]),"cells across",length(unique(X$frame)),"images"))
frames=unique(X$frame)
## Basic plots: frames (i.e. timepoints) per lineage
fpl=grpstats(as.matrix(X$frame), X$lineageId,"numel")
## How many cells from each lineage and track do we have
plyr::count(X$lineageId)
plyr::count(X$trackId)
## How many cells have an assigned parent
unique(X$parentTrackId)
## movie of cells over multiple timepoints
X=X[X$lineageId>=0,] ##Exclude FP detections
plot(X$trackId,X$lineageId)
lims=as.data.frame(apply(X,2,quantile,c(0,1),na.rm=T))
loi=sample(unique(X$lineageId),10); #lineages of interest
for(t in unique(X$frame)){
  X_=X[X$frame==t & X$lineageId %in% loi,]
  plot(X_$Center_of_the_object_0,X_$Center_of_the_object_1,col=X_$lineageId,main=paste("time=",t),xlim=lims$Center_of_the_object_0,ylim=lims$Center_of_the_object_1,pch=20,cex=2)
}


## Add timestamp
FoFs=list.files(A01,pattern=gsub("FoFX","",FoF),full.names = T)
time=sapply(FoFs, function(x) read.table(paste0(x,filesep,"README.md"),sep="\t"))
time=as.data.frame(do.call(rbind,time))
rownames(time)=sapply(rownames(time),function(x) fileparts(x)$name)
time[,2]=sapply(strsplit(time[,2]," "),function(x) paste(x[2:3],collapse = " "))
time[,2]=as.POSIXct(time[,2])
X$POSIXct<-X$FoF<-NA
for(t in 1:nrow(time)){
  ii=which(X$frame==frames[t])
  X$POSIXct[ii]=as.numeric(as.POSIXlt.POSIXct(time[t,2]))
  X$FoF[ii]=rownames(time)[t]
}
X$POSIXct=(X$POSIXct-min(X$POSIXct))/60^2 ## COnvert to hours
plot(X$frame,X$POSIXct)


X=TemporalQuant_sinceANDuntil_division(X)
## Attempt at assigning parent track ID: (TODO -- should be performed by Ilastik's division detection)
X=X[X$lineageId>=0,] 
## For each timepoint except first: 
times=unique(X$POSIXct)
for(it in 2:length(times)){
  t=times[it]
  ## find new cells, C, emerging at that timepoint
  for(i in which(X$first_tpc==t & X$POSIXct==t)){
    X_=X[i,]
    ## find the closest cell, P, from previous timepoint t-1
    d=flexclust::dist2(X_[,XYZ],X[X$POSIXct==times[it-1],XYZ])
    parent=colnames(d)[which.min(d)]
    parent=X[parent,"trackId"]
    ## Assign new track ID to P at timepoint t and forward
    ii = which(X$trackId==parent & X$POSIXct>=t)
    X[ii,"trackId"] = max(X$trackId)+1
    ## Assign P as parent of C and of itself at timepoint t forward
    X$parentTrackId[ii]=parent
    ij = which(X$trackId==X_$trackId & X$POSIXct>=t)
    X$parentTrackId[ij]=parent
  }
}
X=TemporalQuant_sinceANDuntil_division(X)

## Plot division temporal orders
## Classify cell cycle state based on time since last division
X$time_since_division <- X$time_until_division <- NA
# 1) for cells that do have a parent: how long have they been around
ii = which(X$parentTrackId>0)
parents = unique(X$parentTrackId[ii])
X$time_since_division[ii] = X$POSIXct[ii] - X$first_tpc[ii]
# 2) for cells that are a parent: when did they become one?
ii = which(X$trackId %in% parents & X$time_of_division>X$POSIXct)
X$time_until_division[ii] = (X$time_of_division-X$POSIXct)[ii]
ii=which(is.na(X$time_since_division) & !is.na(X$time_until_division))
X$time_since_division[ii]=1.25*max(X$time_since_division,na.rm=T)
## Plots:
plot(X$time_since_division,X$time_until_division)
plot(X$time_since_division, col=X$trackId)
plot(X$time_until_division, col=X$trackId)
## Visualize dividing cells
h5=rhdf5::h5read("./G06_IlastikInput/FoFX007_220523_brightfield.h5",name = FoF)
xshift=150:nrow(h5)
yshift=100:(ncol(h5)-60)
h5=h5[xshift,yshift,,]
coi=X$parentTrackId>0 | X$trackId %in% X$parentTrackId
par(mfrow=c(1,4),mai=c(1,0.25,1,0.1))
for(t in unique(X$frame)[1:4]){
  X_=X[coi & X$frame==t ,]
  image(x=1:nrow(h5),y=1:ncol(h5),z=h5[,,X_$Center_of_the_object_2[1],t+1],col=gray.colors(50),main=paste("time=",t))
  points(X_$Bounding_Box_Maximum_0-xshift[1],X_$Bounding_Box_Minimum_1-yshift[1],col=12+X_$trackId,pch=20,cex=2.5);#,pch=20,xlim=lims$Center_of_the_object_0,ylim=lims$Center_of_the_object_1)
}



## Save time since division info
## Images have been reduced to speedup Ilastik: Scale back up
X_=reverseResize4Ilastik(X)
colnames(X_)[match(XYZ,colnames(X_))]=c("x","y","z")
## Write:
write.table(X_,file=paste0(FROMILASTIK,filesep,FoF,"_DeltaDivision.txt"),sep="\t",row.names = F,quote = F)



## TODO: correlate time_since_division, time_until_division vs. cell size, area, skeweness, etc




## TODO: quantify assymetry in DNA content among daughter cells vs. time until they die or divide
















# 
# X_$group <- balanced_clustering(
#   X_[,c("Center_of_the_object_0","Center_of_the_object_1")],
#   K = round(nrow(X_) / 3) # 5 plants per group
# )