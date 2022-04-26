# source("R/generateImageMask.R")
# FoF="FoF1001_220407_brightfield";
# generateImageMask(FoF)
generateImageMask <- function(FoF, root="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87"){
  library(matlab)
  A05="A05_PostProcessCellposeOutput"
  B06=paste0(root,filesep,"B06_OrganelleMasks")
  success=try(setwd(paste0(root,filesep,A05,"/",FoF,"/All_Cells_coordinates")),silent = T)
  if(class(success)=="try-error"){
    print(paste("No output found for",FoF,"under",A05,". Postprocess segmentation first."))
    return()
  }
  dir.create(paste0(B06,filesep,FoF))
  
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
    iout=paste0(B06,filesep,FoF,filesep,FoF,"_z",z,".tif")
    tiff(iout)
    par(mai=c(0,0,0,0)); image(img,frame.plot=F,axes=F)
    dev.off()
    ## save for gif
    images[[z]]=try(magick::image_read(iout),silent = T)
  }
  
  ## animate at 2 frames per second
  if(all(sapply(images,class)!="try-error")){
    img_joined <- image_join(images)
    img_animated <- image_animate(img_joined, fps = 1)
    image_write(image = img_animated, path = paste0(B06,filesep,FoF,".gif"))
  }else{
    print("Magick not installed. No .gif output generated")
  }
}
