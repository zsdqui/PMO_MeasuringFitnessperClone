Plot_ConcaveHull <- function(xx, yy, zz, lcolor="black", alpha=0.4, add=T, level=0.5/length(xx)) {
  library(MASS)
  ##Remove outliers
  hQ=0.975; lQ=0.025
  iK1=which(xx<=quantile(xx,hQ) & xx>=quantile(xx,lQ))
  iK2=which(yy<=quantile(yy,hQ) & yy>=quantile(yy,lQ))
  iK3=which(zz<=quantile(zz,hQ) & zz>=quantile(zz,lQ))
  iK=intersect(iK1,iK2)
  iK=intersect(iK,iK3)
  xx=xx[iK]; yy=yy[iK]; zz = zz[iK]
  
  ##Contour
  dens2 <- kde3d(xx, yy, zz, lims=c(min(xx)-sd(xx), max(xx)+sd(xx),
                                    min(yy)-sd(yy), max(yy)+sd(yy),
                                    min(zz)-sd(zz), max(zz)+sd(zz) ),n=55  )
  misc3d::contour3d(dens2$d, level=level, dens2$x, dens2$y, dens2$z, color=lcolor, add=add, alpha=alpha); #,drawlabels=F,lwd=2
  # return(cbind(dens2$x,dens2$y, dens2$z))
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