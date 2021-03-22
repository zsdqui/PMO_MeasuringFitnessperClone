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
r3dDefaults$windowRect=c(0,50, 800, 800) 
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getGenesInvolvedIn.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getAllPathways.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
# setwd("~/Projects/PMO/MeasuringFitnessPerClone/code/SingleCellSequencing")
setwd("/mnt/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/SingleCellSequencing")
source("get_compartment_coordinates.R")
source("get_compartment_coordinates_FromAllen.R")
source("alignPathways2Compartments.R")
# source("~/get_compartment_coordinates.R")
# pathwayMapFile = "~/NCBI2Reactome_PE_All_Levels_sapiens.txt"
pathwayMapFile = "NCBI2Reactome_PE_All_Levels_sapiens.txt"
CELLLINE="SNU-668"
ROOTD = tools::file_path_as_absolute("../../data/")
OUTD=paste0("../../results/pathwayCoordinates_3D", filesep, CELLLINE)
dir.create(OUTD,recursive = T)

# coord = get_Compartment_coordinates(300)
coord = get_compartment_coordinates_FromAllen(cytosolF=NULL, nucleusF = paste0(ROOTD,filesep,"3Dbrightfield/allencell/D03_FijiOutput/DNA_anothercell.csv"), mitoF = paste0(ROOTD,filesep,"3Dbrightfield/allencell/D03_FijiOutput/Mito_anothercell.csv"));
# coord = get_compartment_coordinates_FromAllen(nucleusF = paste0(ROOTD,filesep,"3Dbrightfield/allencell/E03_FijiOutput/SNU16C_cell.csv"), mitoF = NULL);
rgl::movie3d(
  movie="CellCompartmentsIn3D_Placeholder", 
  # rgl::spin3d( axis = c(0, 0, 1), rpm = 2),
  rgl::spin3d( axis = c(1, 1, 1), rpm = 3),
  duration = 20, 
  dir = "~/Downloads/",
  type = "gif", 
  clean = TRUE
)
rgl.close()

path2locmap<-read.table(pathwayMapFile, header = FALSE, sep = "\t", dec = ".", comment.char="", quote="", check.names = F, stringsAsFactors = F)
newobject1<-sapply(strsplit(path2locmap$V3, "[", fixed=TRUE), function(x) (x)[2])
path2locmap$V3<-newobject1
path2locmap$V3 <- gsub("]", "", path2locmap$V3)

# Here we replace all compartment names in the path2locmap with their counterpart in coords column names
path2locmap$V3 = str_replace_all(path2locmap$V3, "endoplasmic reticulum membrane", "endoplasmic reticulum")
path2locmap$V3 = str_replace_all(path2locmap$V3, "endoplasmic reticulum lumen", "endoplasmic reticulum")
path2locmap$V3 = str_replace_all(path2locmap$V3, "Golgi membrane", "Gogli apparatus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "Golgi-associated vesicle lumen", "Gogli apparatus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "Golgi-associated vesicle membrane", "Gogli apparatus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "Golgi lumen", "Gogli apparatus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "trans-Golgi network membrane", "Gogli apparatus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "nucleoplasm", "nucleus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "nucleolus", "nucleus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "nuclear envelope", "nucleus")
path2locmap$V3 = str_replace_all(path2locmap$V3, "mitochondrial inner membrane", "mitochondrion")
path2locmap$V3 = str_replace_all(path2locmap$V3, "mitochondrial outer membrane", "mitochondrion")
path2locmap$V3 = str_replace_all(path2locmap$V3, "mitochondrial intermembrane space", "mitochondrion")
path2locmap$V3 = str_replace_all(path2locmap$V3, "mitochondrial matrix", "mitochondrion")
path2locmap$V3 = str_replace_all(path2locmap$V3, "endosome lumen", "endosome")
path2locmap$V3 = str_replace_all(path2locmap$V3, "endosome membrane", "endosome")
path2locmap$V3 = str_replace_all(path2locmap$V3, "late endosome lumen", "endosome")
path2locmap$V3 = str_replace_all(path2locmap$V3, "late endosome membrane", "endosome")
path2locmap$V3 = str_replace_all(path2locmap$V3, "lysomal lumen", "lysosome")
path2locmap$V3 = str_replace_all(path2locmap$V3, "lysomal membrane", "lysosome")
path2locmap$V3 = str_replace_all(path2locmap$V3, "peroxisomal matrix", "peroxisome")
path2locmap$V3 = str_replace_all(path2locmap$V3, "peroxisomal membrane", "peroxisome")

# Expression profiles of all detected genes for clone 9 in the cell line.
cID = 9
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


## Exclude pathways that are active in undefined locations for now. @TODO: map all pathways later
path2locmap = path2locmap[path2locmap$V3 %in% colnames(coord)[apply(coord!=0,2,any)],]
## Exclude pathways that are in endosome or peroxisome: we did not set their coordinates. @TODO later
path2locmap = path2locmap[!path2locmap$V3 %in% c("endosome" ,  "peroxisome"  ),]
## Exclude pathways that are not expressed:
path2locmap = path2locmap[path2locmap$V6 %in% rownames(pq),]
## Rename columns for easier readability
colnames(path2locmap)[c(3,6)]=c("Location","pathwayname")

## Calculate 3D pathway activity maps
LOI=c("nucleus","mitochondrion")
for (cellName in colnames(pq)[1:10]){
  dir.create(paste0(OUTD,filesep,cellName))
  cell1 <- pq[,cellName]
  ## All this info comes from REACTOME + scRNA-seq data:
  pathwayExpressionPerCell = cell1 #this is new
  names(pathwayExpressionPerCell) <- rownames(pq) #this is new
  
  ## The first pathway/location pair we're looking at
  for(j in unique(path2locmap$pathwayname[path2locmap$Location %in% LOI])){
    pmap = cbind(coord, matrix(0,nrow(coord),1))
    colnames(pmap)[ncol(pmap)]=j
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
}

## Animation 3D pathway activity maps
for (cellName in colnames(pq)[1:10]){
  for(outTable in list.files(paste0(OUTD,filesep,cellName), pattern = ".txt", full.names = T )){
    outImage = gsub(".txt",".gif",outTable)
    pmap_ = read.table(file = outTable,sep="\t", header = T, check.names = F, stringsAsFactors = F)
    
    material3d(alpha = 0.1)
    rgl::plot3d(x=pmap_$x, y=pmap_$y, z=pmap_$z,add=F, size=4.91, col=pmap_[,ncol(pmap_)], xlim=quantile(coord$x,c(0,1)), ylim=quantile(coord$y,c(0,1)), zlim=quantile(coord$z,c(0,1)), axes=F, xlab="",ylab="", zlab="", alpha=0.2)
    rgl::movie3d(
      movie=fileparts(outImage)$name, 
      rgl::spin3d( axis = c(1, 1, 1), rpm = 12),
      duration = 4, 
      dir = fileparts(outImage)$path,
      type = "gif", 
      clean = TRUE
    )
    rgl.close()
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


# align VAE output inside compartment
coord_nucl = coord[coord$nucleus==1,]
## Clone_0.0027047_ID107807
vae = read.csv(paste0(ROOTD,"RNAsequencing/A02_210128_VAEoutput/identities_latentSpace3D.csv"))
vae=vae[vae$Location=="nucleus",]
vae=vae[!duplicated(vae),]
scatterplot3d::scatterplot3d(vae$x, vae$y, vae$z, pch=20)
alignPathways2Compartments(coord_nucl[,c("x","y","z")], vae[,c("x","y","z")])


##########################################################################
### Do RNA-seq distributions recapitulate those derived from imaging ? ###
##########################################################################
gatherNucleusCytoExpression <- function(x, path2locmap, q1=0.005, q2=0.995){
  x=x+min(x);
  path2locmap$Location2 = path2locmap$Location
  path2locmap$Location2[path2locmap$Location %in% c("endoplasmic reticulum","mitochondrion" ) ] = "cytosol" 
  
  init = as.matrix(x[1,,drop=F])
  cum = list(cytosol = init, nucleus = init)
  for(p in intersect(rownames(x), path2locmap$pathwayname) ){
    ii = which(path2locmap$pathwayname == p);
    loi = intersect(path2locmap$Location2[ii], names(cum))
    if(!isempty(loi)){
      for(l in loi){
        cum[[l]] = cum[[l]]+x[p,]
      }
    }
  }
  # for(l in names(cum)){
  #   cum[[l]] = cum[[l]] - quantile(cum[[l]],q1)
  #   cum[[l]] = cum[[l]]/quantile(cum[[l]],q2)
  # }
  return(cum)
}


hin = list()
## Imaging data
f=c("~/QuPath/output/NUGC-4A_A8_seedTPd6.txt", "~/QuPath/output/SNU-16C_A7_seedTPd4.txt")
imgdat = lapply(f, function(x) read.table(x,sep="\t", check.names = F, stringsAsFactors = F, header=T));
# names(imgdat) = sapply(f, function(x) fileparts(x)$name)
names(imgdat) = c("NUGC-4","SNU-16")
sapply(imgdat, dim)
for(y in names(imgdat)){
  x = imgdat[[y]]
  x = x[x$`Cell: Area`>=40 & x$`Cell: Area`<=600,]
  x = x[x$`Nucleus: Area`>=10,]
  imgdat[[y]] = x
}

what = "Nucleus/Cell area ratio"
# what = "Nucleus: Area"
# what = "Cell: Area"
# what = "Nucleus: Perimeter"
# what = "Cell: Perimeter"
par(mfrow=c(2,1))
sapply(names(imgdat), function(x) hist(imgdat[[x]][,what], 100, col="blue", border="white", main=x, xlab=what))
hin[["I"]] = as.data.frame(do.call(rbind, lapply(imgdat, function(x) x[,what,drop=F])))
group = sapply(names(imgdat), function(x) rep(x, nrow(imgdat[[x]])))
hin[["I"]]$CellLine = c(group[[1]], group[[2]])
colnames(hin[["I"]])[1] = "Nucleus_Cell_AreaRatio" 
hin[["I"]] = hin[["I"]][hin[["I"]][,1]>= quantile(hin[["I"]][,1],0.00),,drop=F]
hin[["I"]] = hin[["I"]][hin[["I"]][,1]<= quantile(hin[["I"]][,1],1),,drop=F]
easyGgplot2::ggplot2.histogram(data=hin[["I"]], xName = colnames(hin[["I"]])[1], groupName ="CellLine",alpha=0.5, addDensity=T,meanLineColor="white", bins=35, groupColors  = gray.colors(3)[1:2]) +    theme_minimal()

## Sequencing data
load("~/Projects/discussion/grantsAndAwards/Rita Allen 2019/Proposal/SelectionForces_GastricCLs/Results/listOfSeurats_G0G1_p.RData")
seqdat = list(listOfSeurats$`NUGC-4`@assays$RNA@data, listOfSeurats$`KATOIII`@assays$RNA@data)
names(seqdat) = fliplr(names(imgdat))

cytnuc = lapply(seqdat, gatherNucleusCytoExpression, path2locmap, q1=0, q2=1)
names(cytnuc) = names(seqdat)
cytnuc = sapply(cytnuc, function(x) rbind(x$nucleus,x$cytosol))
par(mfrow=c(2,1))
sapply(names(cytnuc), function(x) hist(cytnuc[[x]][1,]/cytnuc[[x]][2,], 100, col="blue", border="white", main=x, xlim=c(0.55,0.75)))

hin[["S"]] = as.data.frame(t(do.call(cbind,cytnuc)))
group = sapply(names(cytnuc), function(x) rep(x, ncol(cytnuc[[x]])))
hin[["S"]]$CellLine = c(group[[1]], group[[2]])
colnames(hin[["S"]])[1:2] = c("Nucleus: Area","Cell: Area" )
# colnames(hin[["S"]])[1:2] = c("Nucleus: Perimeter", "Cell: Perimeter")
hin[["S"]]$Nucleus_Cytosol_ExpressionRatio = hin[["S"]]$`Nucleus: Area`/hin[["S"]]$`Cell: Area`
easyGgplot2::ggplot2.histogram(data=hin[["S"]], xName = "Nucleus_Cytosol_ExpressionRatio", groupName ="CellLine",alpha=0.5, addDensity=T,meanLineColor="white", bins=35, groupColors  = gray.colors(3)[1:2]) +    theme_minimal()

## Merge
tmp = hin$S[,c("Nucleus_Cytosol_ExpressionRatio","CellLine")]
tmp$Assay = "Sequencing"
tmp2 = hin$I[,c("Nucleus_Cell_AreaRatio","CellLine")]
tmp2$Nucleus_Cell_AreaRatio = log(tmp2$Nucleus_Cell_AreaRatio)
tmp2$Nucleus_Cell_AreaRatio = tmp2$Nucleus_Cell_AreaRatio -min(tmp2$Nucleus_Cell_AreaRatio, na.rm=T)
tmp2$Assay = "Imaging"
tmp[,1] = tmp[,1]/mean(tmp[,1])
tmp2[,1] = tmp2[,1]/mean(tmp2[,1])
colnames(tmp) <- colnames(tmp2) <- c("Nucleus_Cytosol", "CellLine", "Assay") 
tmp = rbind(tmp, tmp2)
tmp = tmp[sort(tmp$Assay, index.return=T)$ix,]
par(mai=c(1,3,1,1))
boxplot(Nucleus_Cytosol ~ Assay + CellLine, data=tmp, col=gray.colors(2), horizontal=T, las=2, ylab="")
legend("topleft", c("ScRNA-Seq", "Imaging"), fill=gray.colors(2))
print(plyr::count(paste(tmp$CellLine,tmp$Assay)))


## Co-cluster
cytnucimg <- cytnucseq <- list()
for(what in colnames(hin[["S"]])[1:2] ){
  cytnucseq[[what]] = hin[["S"]][hin[["S"]]$CellLine=="SNU-16", what]
  cytnucseq[[what]] = cytnucseq[[what]]- min(cytnucseq[[what]])
  cytnucseq[[what]] = cytnucseq[[what]]/max(quantile(cytnucseq[[what]],0.995))
  cytnucimg[[what]] = log(imgdat$`SNU-16`[,what])
  cytnucimg[[what]] = cytnucimg[[what]]- quantile(cytnucimg[[what]],0.05)
  cytnucimg[[what]] = cytnucimg[[what]]/max(quantile(cytnucimg[[what]],0.995))
}
par(mfrow=c(2,2))
sapply(names(cytnucimg), function(x) hist(cytnucimg[[x]], xlab=x, main="img",xlim=c(0,1)))
sapply(names(cytnucseq), function(x) hist(cytnucseq[[x]], xlab=x, main="seq",xlim=c(0,1)))

## Sample cells from both assays
N=50
mm = rbind(cbind(cytnucseq[[1]],cytnucseq[[2]]),cbind(cytnucimg[[1]],cytnucimg[[2]]))
rownames(mm) = c(rep("seq",length(cytnucseq[[1]])), rep("img",length(cytnucimg[[1]])))
colnames(mm) = names(cytnucseq)
ii = c(sample(grep("seq",rownames(mm)), N), sample(grep("img",rownames(mm)), N))
mm = as.data.frame(mm[ii,])
mm$ratio=mm[,1]/mm[,2]
mm = mm[is.finite(mm$ratio),]
## Tree
dd = dist(mm[,2:3])
tr = ape::nj(dd)
col = rep("red", length(tr$tip.label))
col[grep("img",tr$tip.label) ] = "blue"
par(mfrow=c(1,1))
plot(tr,show.tip.label = T, tip.color = col, cex=0.6)

## How many sequenced cells match closer to an imaged cell than to a sequenced one:
dd = as.matrix(dd)
dd[dd==0]=NA
rownames(dd) = sapply(strsplit(rownames(dd),".", fixed=T),"[[",1)
colnames(dd) = sapply(strsplit(colnames(dd),".", fixed=T),"[[",1)
bestmatch = paste(rownames(dd), colnames(dd)[apply(dd,1,which.min)])
fr = plyr::count(bestmatch)
rownames(fr) = fr$x
print((fr["img seq",]$freq+fr["seq img",]$freq)/sum(fr$freq))





