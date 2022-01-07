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
source("get_compartment_coordinates_FromAllen.R")
source("alignPathways2Compartments.R")
# source("~/get_compartment_coordinates.R")
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
overlayHist<-function(this,ontothat){
  dat=as.data.frame(c(this,ontothat))
  colnames(dat)="data"
  dat$cat="ontothat"
  dat$cat[1:length(this)]="this"
  p=ggplot(dat,aes(x=data,fill=cat)) + 
    geom_histogram(data=dat, alpha = 0.4)
  return(list(dat=dat,p=p))
}
# this=imgStats[,1]; ontothat=sample(imgStats[,2],100)
overlayDistributions1D<-function(this, ontothat,q=c(0.05,0.5,0.95)){
  # dat=overlayHist(this,ontothat)
  qstat<-function(this, ontothat){
    q=sapply(list(this,ontothat), quantile, q)
    colnames(q)=c("this","ontothat")
    q=as.data.frame(q)
    return(q)
  }
  q1=qstat(this, ontothat)
  ## Shift both medians to 0
  this=this-q1$this[1]
  ontothat=ontothat-q1$ontothat[1]
  ## Overlay upper quantile
  q2=qstat(this, ontothat)
  this=this/(q2$this[3]/q2$ontothat[3])
  
  ## ontothat: Median back to orig
  ontothat=ontothat+q1$ontothat[2]
  this=this+q1$ontothat[2]
  dat=overlayHist(this,ontothat)
  return(dat)
}
# pathwayMapFile = "~/NCBI2Reactome_PE_All_Levels_sapiens.txt"
pathwayMapFile = "NCBI2Reactome_PE_All_Levels_sapiens.txt"
CELLLINE="NCI-N87"
ROOTD = tools::file_path_as_absolute("../../data/")
OUTD=paste0("../../results/pathwayCoordinates_3D", filesep, CELLLINE)
dir.create(OUTD,recursive = T)


############################################################
### Load coordinates of various compartments for one cell ##
# coord = get_Compartment_coordinates(300)
# coord = get_compartment_coordinates_FromAllen(cytosolF=NULL, nucleusF = paste0(ROOTD,filesep,"3Dbrightfield/allencell/D03_FijiOutput/DNA_anothercell.csv"), mitoF = paste0(ROOTD,filesep,"3Dbrightfield/allencell/D03_FijiOutput/Mito_anothercell.csv"));
coord = get_compartment_coordinates_FromAllen(nucleusF = "../../data/3Dbrightfield/allencell/G03_CellposeOutput/0_prediction_c0.model.p_cell_109_coordinates.csv", mitoF = "../../data/3Dbrightfield/allencell/G03_CellposeOutput/0_prediction_c0.model.p_cell_108_coordinates.csv", XYZCOLS = c("x","y","z"));
rgl::movie3d(
  movie="CellCompartmentsIn3D_Placeholder", 
  rgl::spin3d( axis = c(1, 1, 1), rpm = 3),
  duration = 20, 
  dir = "~/Downloads/",
  type = "gif", 
  clean = TRUE
)
rgl.close()
## Calculate volume
hull <- convhulln(coord[coord$nucleus==1,1:3], options = "FA")
print(hull$vol)

# Expression profiles of all detected genes for clone 9 in the cell line.
cID = 3
clones = cloneid::getSubclones(CELLLINE, whichP="TranscriptomePerspective")
pqFile = paste0(OUTD,filesep,names(clones)[cID],".RObj")
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
  save(pq, file=paste0(OUTD,filesep,names(clones)[cID],".RObj"))
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
par(mfrow=c(2,2))
sapply(colnames(seqStats), function(x) hist(seqStats[,x],xlab=x))
## read stats from 3D imaging and plot side by side:
imgStats=read.table("../../data/3Dbrightfield/allencell/G04_segmentationStats/0_prediction_c0.model_stats.txt",sep="\t",check.names = F,stringsAsFactors = F)
sapply(colnames(imgStats), function(x) hist(imgStats[,x],xlab=x))
## Overlay distributions
seqStatsMapped=list()
for(i in 1:ncol(imgStats)){
  overlayHist(seqStats[,i],imgStats[,i])$p
  seqStatsMapped[[colnames(imgStats)[i]]]=overlayDistributions1D(seqStats[,i],imgStats[,i], q=c(0.1,0.5,0.9))
}
seqStatsMapped$volume$p 
ggsave(filename = "~/Downloads/volumeFeature.png",width = 4,height = 3)
## Overwrite seqStats
seqStats=as.data.frame(sapply(seqStatsMapped,function(x) x$dat$data))
colnames(seqStats)=colnames(imgStats)
seqStats=seqStats[sample(nrow(seqStats),size = nrow(imgStats)),]
## co-cluster image and sequencing stats
imgStats$type="img"
seqStats$type="seq"
##@TODO: before reassigning colnames, first ensure they are in desired order (which pair of features match)
stats=rbind(imgStats,seqStats)
rownames(stats) = paste(rownames(stats), stats$type)
dd = dist(stats[,1:2])
tr = ape::nj(dd)
plot(tr)
col = rep("red", length(tr$tip.label))
col[grep("img",tr$tip.label) ] = "blue"
par(mfrow=c(1,1))
plot(tr,show.tip.label = T, tip.color = col, cex=0.36)
legend("topright",c("sequencing","imaging"),fill=c("red","blue"))


## locations per pathway:
lpp = sapply(unique(path2locmap$pathwayname), function(x) unique(path2locmap$Location[path2locmap$pathwayname==x]))
# lpp = lpp[sample(length(lpp),10)]; ## use only subset for testing
lpp = lpp[grep("Cycle",names(lpp),value=T)[c(7:9)]]
save(file='~/Downloads/tmp_coord.RObj', list=c('coord','OUTD', 'lpp','pq','path2locmap'))

## Calculate 3D pathway activity maps
pathwayColors=rainbow(length(lpp))
names(pathwayColors)=names(lpp)
LOI=c("nucleus","mitochondrion")
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
    # In this step we populate our pmap with our randomly selected x
    # coordinate and it's matching y coordiate from the coord object
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
  for(compartment in LOI){
    ii=which(coord[,compartment]==1)
    hull=Plot_ConcaveHull(coord[ii,1], coord[ii,2], coord[ii,3], lcolor ="gray", alpha=0.075)
  }
  for(j in names(lpp)){
    pmap_ = pmap[pmap[,j]>0,]
    rgl::points3d(x=pmap_$x, y=pmap_$y, z=pmap_$z,add=F, size=6,col=pathwayColors[j], xlim=quantile(coord$x,c(0,1)), ylim=quantile(coord$y,c(0,1)), zlim=quantile(coord$z,c(0,1)), axes=F, xlab="",ylab="", zlab="", alpha=0.3)
  }
  legend3d("topright",names(pathwayColors),fill=pathwayColors)
  
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

