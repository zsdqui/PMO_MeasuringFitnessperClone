color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut))/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut))) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


Plot_ConcaveHull <- function(xx, yy, zz, lcolor="black", alpha=0.4, add=T, level=0.5/length(xx),xlim=quantile(xx,c(0,1)),ylim=quantile(yy,c(0,1)),zlim=quantile(zz,c(0,1)), other=NULL, png=NULL) {
  library(MASS)
  # print("datapoints nucleus:"); print(length(xx))
  # print("datapoints mito:"); print(nrow(other))
  
  ##Remove outliers
  # hQ=0.975; lQ=0.025
  hQ=1; lQ=0
  iK1=which(xx<=quantile(xx,hQ) & xx>=quantile(xx,lQ))
  iK2=which(yy<=quantile(yy,hQ) & yy>=quantile(yy,lQ))
  iK3=which(zz<=quantile(zz,hQ) & zz>=quantile(zz,lQ))
  iK=intersect(iK1,iK2)
  iK=intersect(iK,iK3)
  xx=xx[iK]; yy=yy[iK]; zz = zz[iK]
  qx=quantile(xx,c(0,1))
  qy=quantile(yy,c(0,1))
  qz=quantile(zz,c(0,1))
  
  ##Contour
  n=1+c(qx[2]-qx[1], qy[2]-qy[1], qz[2]-qz[1])
  # print("nucleus:"); print(n)
  dens2 <- misc3d::kde3d(xx, yy, zz,n=n ) 
  RES=c(x=dens2$x[2]-dens2$x[1],y=dens2$y[2]-dens2$y[1],z=dens2$z[2]-dens2$z[1])
  
  ## Embed in larger array to render plots comparable across cells
  dens <- array(0,dim=ceil(cbind(xlim,ylim,zlim)[2,]/RES))
  # print("all:"); 
  print(dim(dens))
  dens[ceil(dens2$x/RES["x"]),ceil(dens2$y/RES["y"]), ceil(dens2$z/RES["z"])]=dens2$d
  ## Plot
  
  if(!is.null(png)){
    png(png)
  }
  plot.new()
  par(mai=c(0,0,0,0))
  plot3D::isosurf3D(x=1:nrow(dens),y=1:ncol(dens),z=1:dim(dens)[3],colvar = dens, col=lcolor, add=add, alpha=alpha,bty="n",phi = 40, theta = 40); #,drawlabels=F,lwd=2
  # misc3d::contour3d(dens, level=level, color=lcolor, add=add, alpha=alpha); #,drawlabels=F,lwd=2
  
  ## add other compartment to coordinate system:
  dens_o=NULL
  if(!is.null(other)){
    other=sweep(other,MARGIN = 2, STATS = RES, FUN = "/")
    plot3D::scatter3D(other[,1], other[,2], other[,3],pch3d=20, size=5, axes=F, xlab="",ylab="", zlab="",col="red",alpha=0.04, add=T,xlim=xlim,ylim=ylim,zlim=zlim)
    dens_o=dens
    dens_o[T]=0
    for(i in 1:nrow(other)){
      dens_o[other[i,1], other[i,2], other[i,3]] =1
      ##@TODO: should this be assigned a smaller value?
    }
  }
  if(!is.null(png)){
    dev.off()
  }
  # # ##TESTPLOT
  # xyz=varbvs::grid3d(1:dim(dens)[1],1:dim(dens)[2],1:dim(dens)[3]);
  # xyz=do.call(cbind,xyz)
  # out=cbind(xyz,hull[T]);
  # out=out[out[,4]>1E-7,]
  # plot.new()
  # plot3D::scatter3D(out[,1], out[,2], out[,3],pch=20, size=0.2, axes=F, xlab="",ylab="", zlab="",alpha=0.04, add=F,colvar =out[,4],bty="n", colkey=F, type = "p")
  
  
  return(list(nucleus=dens, other=dens_o))
}


get_compartment_coordinates_FromAllen <-function(cytosolF="~/Downloads/fijitestout2.csv", mitoF="~/Downloads/fijitestout2.csv", nucleusF="~/Downloads/fijitestout3.csv",XYZCOLS = c("CX..pix.", "CY..pix.", "CZ..pix."), size=0.01){
  COMPCOLS = c("cytosol", "endoplasmic reticulum", "Golgi apparatus", "nucleus", "mitochondrion", "endosome", "lysosome", "peroxisome")
  alpha=rep(1,length(COMPCOLS));
  names(alpha)=COMPCOLS
  alpha["nucleus"]=1
  
  dummy =as.data.frame( matrix(0,1,length(COMPCOLS)))
  colnames(dummy) = COMPCOLS
  
  dummy[rep(seq_len(nrow(dummy)), each = 2), ]
  
  nucl=read.csv(nucleusF)
  nucl = cbind(nucl[,XYZCOLS], dummy[rep(seq_len(nrow(dummy)), each = nrow(nucl)), ])
  nucl$nucleus = 1;
  coord = nucl
  if(!is.null(mitoF)){
    mito=read.csv(mitoF)
    mito = cbind(mito[,XYZCOLS], dummy[rep(seq_len(nrow(dummy)), each = nrow(mito)), ])
    mito$`mitochondrion`=1
    coord = rbind(coord, mito);
  }
  if(!is.null(cytosolF)){
    cyto=read.csv(cytosolF)
    cyto = cbind(cyto[,XYZCOLS], dummy[rep(seq_len(nrow(dummy)), each = nrow(cyto)), ])
    cyto$cytosol = 1;
    coord = rbind(coord, cyto)
  }
  colnames(coord)[1:3] = c("x","y","z")
  
  ## draw cell to check coordinate assignment
  library(RColorBrewer)
  colormap = brewer.pal(length(COMPCOLS),"Paired")
  names(colormap) = COMPCOLS
  tmp = quantile(coord$z,c(0,1))
  space=tmp[2]-tmp[1]
  scatterplot3d::scatterplot3d(coord$x, coord$y, coord$z, pch=20,cex.symbols = 0.1, zlim=c(tmp[1]-space/1.5, tmp[2]+space/1.5))
  rgl::plot3d(coord$x, coord$y, coord$z, pch=20, zlim=c(tmp[1]-space/2, tmp[2]+space/2), size=size, axes=F, xlab="",ylab="", zlab="")
  for(o in names(colormap)){
    print(o)
    coord_ = coord[coord[,o]==1,]
    rgl::points3d(coord_$x, coord_$y, coord_$z,col=colormap[o],add=T, size=size,alpha=alpha[o])
    print(nrow(coord_))
  }
  rgl::legend3d("topleft", names(colormap), fill=colormap, bty='n',cex=1.7)
  
  return(coord)
}

dirCreate<-function(dname, recursive=T, permission=NULL){
  dir.create(dname, recursive = recursive)
  if(!is.null(permission)){
    system(paste("chmod -R", permission, dname))
  }
}



resize4Ilastik<-function(img, xydim = 255){
  img=EBImage::resize(img,h = xydim, w=xydim)
  return(img)
}


##@TODO: test
reverseResize4Ilastik<-function(df, xydim_from=255, xydim_to = 1024){
  XY=c("Center_of_the_object_0","Center_of_the_object_1")
  fac=xydim_to/xydim_from
  df[,XY]=df[,XY]*fac
  return(df)
}


collectTensorsAsVectors<-function(FoF, indir="A06_multiSignals_Linked", root="~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87"){
  f=list.files(paste0(root,filesep,indir,filesep,FoF),pattern=".rds", full.names = T)
  f_m=grep("mito",f,value=T)
  f_n=grep("nuc",f,value=T)
  mat=c()
  for(i in 1:length(f_n)){
    x=readRDS(f_n[i])
    x_m=readRDS(f_m[i])
    x_m[x_m==1]=max(x)*2
    jj=which(x==0)
    x[jj]=x_m[jj]
    mat=rbind(mat,as.numeric(x))
  }
  rownames(mat)=sapply(f_n, function(x) strsplit(fileparts(x)$name,"_")[[1]][2])
  return(mat)
}


getTimeStampsFromMetaData<-function(FoFs, root="A01_rawData", xmlfiles=NULL){
  dat=list()
  for(FoF in FoFs){
    if(is.null(xmlfiles)){
      f=list.files(paste0(root,filesep,FoF,filesep,"MetaData"),pattern="_Properties.xml", full.names = T)
    }else{
      f=grep(FoF,xmlfiles,value = T)
    }
    x=grep('StartTime',readLines(f),value=T)
    x=strsplit(x,">")[[1]][2]
    x=strsplit(x,"<")[[1]][1]
    x=strsplit(x," ")[[1]]
    time=strsplit(x[2],":")[[1]]
    hour=as.numeric(time[1])
    if(x[3]=="PM"){
      hour=hour+12
    }
    time[1]=hour
    time=paste(time,collapse = ":")
    date=strsplit(x[1],"/")[[1]]
    date=paste(date[c(3,1,2)],collapse = "-")
    dat[[FoF]]=as.POSIXct(paste0(date," ",time))
  }
  dat=sapply(dat, function(x) as.numeric(as.POSIXlt.POSIXct(x)))
  dat=(dat-min(dat))/60^2 ## COnvert to hours
  
  return(dat)
}

getZslice<-function(FoF, slice, root="A06_multiSignals_Linked",plot=T, signal="nucleus.p"){
  csv=list.files(paste0(root,filesep,FoF),pattern=signal,full.names = T)
  out=list()
  for (x in csv){
    y=read.csv(x)
    y$cell=fileparts(x)$name
    y$ID=as.numeric(strsplit(fileparts(x)$name,"_")[[1]][3])
    y=y[y$z==slice,,drop=F]
    out[[fileparts(x)$name]]=y
  }
  print(paste(sum(sapply(out,nrow)==0),"cells don't have representation on this z-slice"))
  out=do.call(rbind, out)
  if(plot){
    try(plot(out$x,out$y,col=out$ID))
  }
  return(out)
}
