# source("R/generateImageMask.R")
# FoF="FoF1001_220407_brightfield";
# generateImageMask(FoF)
generateImageMask <- function(FoF, INDIR="A05_PostProcessCellposeOutput", OUTDIR="B06_OrganelleMasks", root="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87"){
  library(matlab)
  OUTDIR=paste0(root,filesep,OUTDIR)
  olddir=getwd()
  success=try(setwd(paste0(root,filesep,INDIR,"/",FoF,"/All_Cells_coordinates")),silent = T)
  if(class(success)=="try-error"){
    print(paste("No output found for",FoF,"under",INDIR,". Postprocess segmentation first."))
    return()
  }
  dir.create(paste0(OUTDIR,filesep,FoF))
  
  tmp=list.files()
  cells=1:length(tmp)
  names(cells)=tmp
  zstack=70
  images=list()
  for(z in 1:zstack){
    img=matrix(NA,1024,1024)
    for(cell in names(cells)){
      csv=read.csv(cell)
      for(i in which(csv$z==z)){
        img[csv$x[i],csv$y[i]]=cells[cell]
      }
    }
    
    ## write image
    z_=z
    if(z<10){
      z_=paste0("0",z)
    }
    iout=paste0(OUTDIR,filesep,FoF,filesep,FoF,"_z",z_,".tif")
    tiff(iout)
    par(mai=c(0,0,0,0)); image(img,frame.plot=F,axes=F)
    dev.off()
    ## save for gif
    images[[z]]=try(magick::image_read(iout),silent = T)
  }
  
  ## animate at 2 frames per second
  cmd=paste0("convert -delay 20 -loop 0 ",OUTDIR,filesep,FoF,filesep,"*.tif ",OUTDIR,filesep,FoF,".gif ")
  system(cmd)
  
  setwd(olddir)
  # if(all(sapply(images,class)!="try-error")){
  #   img_joined <- image_join(images)
  #   img_animated <- image_animate(img_joined, fps = 1)
  #   image_write(image = img_animated, path = paste0(OUTDIR,filesep,FoF,".gif"))
  # }else{
  #   print("Magick not installed. No .gif output generated")
  # }
}
