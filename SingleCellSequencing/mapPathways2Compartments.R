options(java.parameters = "-Xmx7g")
# Load libraries and devtools
library(matlab)
library(GSVA)
library(cloneid)
library(ggplot2)
library(dplyr)
library(plyr)
library(tibble)
library(stringr)
library(scales)
library(rgl)
library(misc3d)
library(geometry)
r3dDefaults$windowRect=c(0,50, 800, 800) 
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getGenesInvolvedIn.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getAllPathways.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
setwd("~/Projects/PMO/MeasuringFitnessPerClone/code/SingleCellSequencing")
# setwd("/mnt/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/SingleCellSequencing")
source("get_compartment_coordinates.R")
source("../3D_Imaging/R/Utils.R")
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
  return(dens2)
}
# pathwayMapFile = "~/NCBI2Reactome_PE_All_Levels_sapiens.txt"
pathwayMapFile = "NCBI2Reactome_PE_All_Levels_sapiens.txt"
CELLLINE="NCI-N87"
ROOTD = tools::file_path_as_absolute("../../data/")
OUTD=paste0("../../results/pathwayCoordinates_3D", filesep, CELLLINE)
B01=paste0("../../data/GastricCancerCL/RNAsequencing/B01_220112_pathwayActivity", filesep, CELLLINE)
B02=paste0("../../data/GastricCancerCL/RNAsequencing/B02_220112_seqStats", filesep, CELLLINE)
dir.create(OUTD,recursive = T)
LOI=c("gray","pink","cyan")
alpha=rep(0.05,length(LOI));
names(alpha) <- names(LOI) <- c("nucleus","mitochondrion","cytosol")
# cellID=251
cellID=c(37,35,15)[1]
FoF="../../data/GastricCancerCL/3Dbrightfield/NCI-N87/A06_multiSignals_Linked/FoF13_220228_fluorescent.cytoplasm/"

############################################################
### Load coordinates of various compartments for one cell ##
# coord = get_Compartment_coordinates(300)
rgl::close3d()
coord = get_compartment_coordinates_FromAllen(nucleusF = paste0(FoF,"nucleus.p_cell_",cellID,"_coordinates.csv"), 
                                              mitoF = paste0(FoF,"mito.p_cell_",cellID,"_coordinates.csv"),
                                              cytosolF = paste0(FoF,"cytoplasm.t_cell_",cellID,"_coordinates.csv"), XYZCOLS = c("x","y","z"), size = 1);
# coord = get_compartment_coordinates_FromAllen(nucleusF = paste0("../../data/GastricCancerCL/3Dbrightfield/NCI-N87/H05_multiOrganelles_Linked/nucleus.p_cell_",cellID,"_coordinates.csv"), mitoF = paste0("../../data/GastricCancerCL/3Dbrightfield/NCI-N87/H05_multiOrganelles_Linked/mito.p_cell_",cellID,"_coordinates.csv"), XYZCOLS = c("x","y","z"), size = 3);
rgl::movie3d(
  movie="~/Downloads/CellCompartmentsIn3D_Placeholder", 
  rgl::spin3d( axis = c(1, 1, 1), rpm = 18),
  duration = 3, 
  dir = "~/Downloads/",
  type = "gif", 
  clean = TRUE
)
rgl.close()
## Calculate volume
hull <- convhulln(coord[coord$nucleus==1,1:3], options = "FA")
print(hull$vol)
## Plot organelle hulls
rgl::close3d()
for(compartment in names(LOI)){
  ii=which(coord[,compartment]==1)
  hull=Plot_ConcaveHull(coord[ii,1], coord[ii,2], coord[ii,3], lcolor =LOI[compartment], alpha=alpha[compartment])
}
legend3d("topleft",names(LOI),fill=LOI)


# Expression profiles of all detected genes for clone 9 in the cell line.
cID = 2
clones = cloneid::getSubclones(CELLLINE, whichP="TranscriptomePerspective")
pqFile = paste0(B01,filesep,names(clones)[cID],".RObj")
if(file.exists(pqFile)){
  load(pqFile)
}else{
  p = getSubProfiles(as.numeric(cloneid::extractID(names(clones)[cID])))
  # Exclude copy number profile (keep only expression profile)
  ex = p[-grep(":", rownames(p)),]
  # Now let's quantify pathway expression based on this expression
  gs=getAllPathways(include_genes=T, loadPresaved = T);     
  gs=gs[sapply(gs, length)>=5]
  pq <- gsva(ex, gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)
  pq <- rescale(pq, to = c(0,30000))
  save(pq, file=paste0(B01,filesep,names(clones)[cID],".RObj"))
}
ccState = sapply(colnames(pq), function(x) cloneid::getState(x,whichP = "TranscriptomePerspective"))


## Pathway information:
path2locmap<-read.table(pathwayMapFile, header = FALSE, sep = "\t", dec = ".", comment.char="", quote="", check.names = F, stringsAsFactors = F)
# Here we replace all compartment names in the path2locmap with their counterpart in coords column names
path2locmap$V3 <- sapply(strsplit(path2locmap$V3, "[", fixed=TRUE), function(x) (x)[2])
path2locmap$V3 <- gsub("]", "", path2locmap$V3)
aliasMap = list("endoplasmic reticulum membrane"="endoplasmic reticulum", "endoplasmic reticulum lumen" = "endoplasmic reticulum", "Golgi membrane" = "Gogli apparatus", "Golgi-associated vesicle lumen" = "Gogli apparatus", "Golgi-associated vesicle membrane" = "Gogli apparatus", "Golgi lumen" = "Gogli apparatus", "trans-Golgi network membrane" = "Gogli apparatus", "nucleoplasm" = "nucleus", "nucleolus" = "nucleus", "nuclear envelope" = "nucleus", "mitochondrial inner membrane" = "mitochondrion", "mitochondrial outer membrane" = "mitochondrion", "mitochondrial intermembrane space" = "mitochondrion", "mitochondrial matrix" = "mitochondrion", "endosome lumen" = "endosome", "endosome membrane" = "endosome", "late endosome lumen" = "endosome", "late endosome membrane" = "endosome", "lysomal lumen" = "lysosome", "lysomal membrane" = "lysosome", "peroxisomal matrix" = "peroxisome", "peroxisomal membrane" = "peroxisome")
for(x in names(aliasMap)){
  path2locmap$V3 = str_replace_all(path2locmap$V3, x, aliasMap[[x]])
}
## Exclude pathways that are active in undefined locations for now. @TODO: map all pathways later
path2locmap = path2locmap[path2locmap$V3 %in% colnames(coord)[apply(coord!=0,2,any)],]
## Exclude pathways that are in endosome or peroxisome: we did not set their coordinates. @TODO later
path2locmap = path2locmap[!path2locmap$V3 %in% c("endosome" ,  "peroxisome"  ),]
## Exclude pathways that are not expressed:
path2locmap = path2locmap[path2locmap$V6 %in% rownames(pq),]
## Rename columns for easier readability
colnames(path2locmap)[c(3,6)]=c("Location","pathwayname")



## calculate pq stats: how much is going on in one compartment vs. another?
seqStats=lapply(unique(path2locmap$Location), function(x) pq[rownames(pq) %in% path2locmap[path2locmap$Location==x,"pathwayname"],])
names(seqStats)=unique(path2locmap$Location)
seqStats=cbind(sapply(seqStats,colMeans),sapply(seqStats,colSums))
colnames(seqStats)=paste(c("meanE","sumE"),colnames(seqStats))
seqStats=as.data.frame(seqStats)
seqStats=cbind(seqStats,t(pq[grep("Cell cycle",rownames(pq),ignore.case = T),]))
write.table(seqStats, paste0(B02,filesep,fileparts(pqFile)$name,".txt"),sep="\t")


## locations per pathway:
lpp = sapply(unique(path2locmap$pathwayname), function(x) unique(path2locmap$Location[path2locmap$pathwayname==x]))
# lpp = lpp[sample(length(lpp),10)]; ## use only subset for testing
lpp = lpp[grep("Cycle",names(lpp),value=T,ignore.case = T)]#[c(7:9)]]
lpp=lpp[order(sapply(lpp,paste,collapse=","),decreasing = F)]
lpp=lpp[1:10]
save(file='~/Downloads/tmp_coord.RObj', list=c('coord','OUTD', 'lpp','pq','path2locmap'))

## Calculate 3D pathway activity maps
pathwayColors=rainbow(length(lpp))
names(pathwayColors)=names(lpp)
alpha["nucleus"]=1
alpha["cytoplasm"]=0.01
for (cellName in colnames(pq)[1]){
  dir.create(paste0(OUTD,filesep,cellName))
  pathwayExpressionPerCell <- pq[,cellName]/20
  names(pathwayExpressionPerCell) <- rownames(pq) 
  pmap = cbind(coord, matrix(0,nrow(coord),length(lpp)))
  colnames(pmap)[(length(coord)+1):ncol(pmap)] = names(lpp)
  
  ## The first pathway/location pair we're looking at
  for(j in names(lpp)){
    outImage = paste0(OUTD,filesep,cellName,filesep,gsub(" ","_",gsub("/","-",j,fixed = T)),".gif")
    outTable = paste0(OUTD,filesep,cellName,filesep,gsub(" ","_",gsub("/","-",j,fixed = T)),".txt")
    if(file.exists(outImage)){
      print(paste("Skipping",j,"because image already saved"))
      next;
    }
    P = path2locmap[path2locmap$pathwayname==j,,drop=F]
    fr =  plyr::count(P$Location)
    fr$freq = fr$freq/sum(fr$freq)
    rownames(fr) = fr$x
    ## Candidates of indices
    idx_Candid = lapply(rownames(fr), function(location) which(coord[,location]==1) )
    names(idx_Candid) = rownames(fr)
    # Here we take a random sample of x coordinates for our pathway 
    idx = lapply(names(idx_Candid), function(x) sample(idx_Candid[[x]], pathwayExpressionPerCell[j]*fr[x,"freq"], replace = T)  )
    names(idx) = names(idx_Candid)
    # And here we tally how many times each x coordinate appears in the sampling
    idx = plyr::count(unlist(idx))
    # In this step we populate our pmap with our randomly selected x coordinate and it's matching y coordiate from the coord object
    for(i in 1:nrow(idx)){
      # pmap[coord$x[idx$x[i]], coord$y[idx$x[i]]] = idx$freq[i]
      pmap[idx$x[i],j] =  idx$freq[i]
    }
    ## Print statement
    print(paste("Processed pathway",j))
    
    # Image of the pathway map 
    pmap_ = pmap[pmap[,j]>0,]
    write.table(pmap_[,c("x","y","z",j)], file = outTable,sep="\t",quote = F, row.names = F)
    # # Image of the pathway map 
    # png(outImage,width = 400, height = 400)
    # image(pmap[j,,],col =rainbow(100),xaxt = "n",yaxt = "n"); #,main=paste(j,cloneid::extractID(cellName))
    # dev.off()
  }
  ##  Plot all pathways together
  rgl::close3d()
  for(compartment in names(LOI)){
    ii=which(coord[,compartment]==1)
    hull=Plot_ConcaveHull(coord[ii,1], coord[ii,2], coord[ii,3], lcolor =LOI[compartment], alpha=alpha[compartment])
  }
  legend3d("topleft",names(LOI),fill=LOI)
  for(j in names(lpp)){
    pmap_ = pmap[pmap[,j]>0,]
    rgl::points3d(x=pmap_$x, y=pmap_$y, z=pmap_$z,add=F, size=5,col=pathwayColors[j], xlim=quantile(coord$x,c(0,1)), ylim=quantile(coord$y,c(0,1)), zlim=quantile(coord$z,c(0,1)), axes=F, xlab="",ylab="", zlab="", alpha=0.7)
  }
  legend3d("topright",names(lpp),fill=pathwayColors[names(lpp)])
  
}

## Animation 3D pathway activity maps
load('~/Downloads/tmp_coord.RObj')
detach('package:GSVA', unload=TRUE)
library(rgl)
library(magick)
for (cellName in list.dirs(OUTD, full.names = F)){
  print(cellName)
  for(outTable in list.files(paste0(OUTD,filesep,cellName), pattern = ".txt", full.names = T )){
    outImage = gsub(".txt",".gif",outTable)
    print(outImage)
    if(!file.exists(outImage)){
      pmap_ = read.table(file = outTable,sep="\t", header = T, check.names = F, stringsAsFactors = F)
      
      r3dDefaults$windowRect = c(0,0,700,700)
      rgl::material3d(alpha = 0.1)
      rgl::points3d(x=pmap_$x, y=pmap_$y, z=pmap_$z,add=F, size=4.91, col=pmap_[,ncol(pmap_)], xlim=quantile(coord$x,c(0,1)), ylim=quantile(coord$y,c(0,1)), zlim=quantile(coord$z,c(0,1)), axes=F, xlab="",ylab="", zlab="", alpha=0.2)
      Sys.sleep(5)
      try(rgl::movie3d(
        movie=matlab::fileparts(outImage)$name, 
        rgl::spin3d( axis = c(1, 1, 1), rpm = 12),
        duration = 8, 
        dir = matlab::fileparts(outImage)$path,
        type = "gif", 
        clean = T
      ))
      rgl::rgl.close()
      
      ##################
      #### clean up ####
      ##################
      detach('package:rgl', unload=TRUE)
      library(crosstalk)
      library(manipulateWidget)
      library(miniUI)
      library(shiny)
      library(shinythemes)
      detach('package:crosstalk', unload=TRUE)
      detach('package:manipulateWidget', unload=TRUE)
      detach('package:miniUI', unload=TRUE)
      detach('package:shinythemes', unload=TRUE)
      detach('package:shiny', unload=TRUE)
      detach("package:magick", unload = TRUE)
      ## Remove cache generated by magick package
      f=list.files('/tmp/Rtmp2quoCU', pattern='magick', full.names=T); 
      if(length(f)>50){
        for (x in f){ 
          file.remove(x)
        }
      }
    }
  }
}




# f=`ls`
# for x in $f; do
#   echo $x;
#   tar -czvf $x.tar.gz $x;
# done


# # Finally we produce an image of the pathway map for testing:
# library(matlab)
# library(cloneid)
# of = list.files(path = "~/Downloads", pattern = "pathwayCoordinatesStack_Clone",full.names = T)[4]
# pmap2D = read.table(file=of, sep="\t", check.names = F, strip.white = F, header = T)
# par(mfrow=c(2,2));
# for(p in c("Mitotic Metaphase and Anaphase","Apoptosis","Metabolism","Mitochondrial protein import")){
#   tmp=pmap2D[pmap2D$pathway==p,]
#   pdf(paste0("~/Downloads/",p,".pdf"))
#   image(as.matrix(tmp[,-1]),main=paste(p,cloneid::extractID(fileparts(of)$name)),col =rainbow(100))
# }
# dev.off()


# # align VAE output inside compartment
# coord_nucl = coord[coord$nucleus==1,]
# ## Clone_0.0027047_ID107807
# vae = read.csv(paste0(ROOTD,"RNAsequencing/A02_210128_VAEoutput/identities_latentSpace3D.csv"))
# vae=vae[vae$Location=="nucleus",]
# vae=vae[!duplicated(vae),]
# scatterplot3d::scatterplot3d(vae$x, vae$y, vae$z, pch=20)
# alignPathways2Compartments(coord_nucl[,c("x","y","z")], vae[,c("x","y","z")])

