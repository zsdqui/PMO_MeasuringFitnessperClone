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

