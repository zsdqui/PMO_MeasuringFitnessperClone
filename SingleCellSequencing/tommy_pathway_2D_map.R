# You must always run this when you first start R or you will not have enough memory available when you need to ping CloneID
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
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getGenesInvolvedIn.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getAllPathways.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
source("~/get_compartment_coordinates.R")
# source("~/Projects/PMO/MeasuringFitnessPerClone/code/get_compartment_coordinates.R")

CELLLINE="SNU-668"

coord = get_Compartment_coordinates()
colnames(coord) = c("x","y","cytosol", "endoplasmic reticulum", "Golgi apparatus", "nucleus", "mitochondrion", "endosome", "lysosome", "peroxisome")

path2locmap<-read.table("~/NCBI2Reactome_PE_All_Levels_sapiens.txt", header = FALSE, sep = "\t", dec = ".", comment.char="", quote="", check.names = F, stringsAsFactors = F)
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

# Expression profiles of all detected genes for clone 7 in the cell line.
clones = cloneid::getSubclones(CELLLINE, whichP="TranscriptomePerspective")
p = getSubProfiles(as.numeric(cloneid::extractID(names(clones)[7])))
# Exclude copy number profile (keep only expression profile)
ex = p[-grep(":", rownames(p)),]
# Now let's quantify pathway expression based on this expression
gs=getAllPathways(include_genes=T, loadPresaved = T);     
gs=gs[sapply(gs, length)>=5]
pq <- gsva(ex, gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)
# Now we're going to look at one cell from all the transcriptomes we grabbed, and un-list it
pq <- rescale(pq, to = c(0,600))
cellName =colnames(pq)[1]
cell1 <- pq[,cellName]

## All this info comes from REACTOME + scRNA-seq data.
pathwayExpressionPerCell = cell1 #this is new
names(pathwayExpressionPerCell) <- rownames(pq) #this is new

## Exclude pathways that are active in undefined locations for now. @TODO: map all pathways later
path2locmap = path2locmap[path2locmap$V3 %in% colnames(coord),]
## Exclude pathways that are in endosome or peroxisome: we did not set their coordinates. @TODO later
path2locmap = path2locmap[!path2locmap$V3 %in% c("endosome" ,  "peroxisome"  ),]
## Exclude pathways that are not expressed:
path2locmap = path2locmap[path2locmap$V6 %in% rownames(pq),]
## Rename columns for easier readability
colnames(path2locmap)[c(3,6)]=c("Location","pathwayname")

pmap = array(0, dim=c(length(rownames(pq)),nrow(coord), nrow(coord)))
rownames(pmap) = rownames(pq)
## The first pathway/location pair we're looking at
for(j in unique(path2locmap$pathwayname)){
  P = path2locmap[path2locmap$pathwayname==j,,drop=F]
  idx_Candid = which(apply(coord[,P$Location,drop=F]==1,1, any)); ## Candidates of indices
  # Here we take a random sample of x coordinates for our pathway 
  idx = sample(idx_Candid, pathwayExpressionPerCell[P$pathwayname], replace = T) 
  # And here we tally how many times each x coordinate appears in the sampling
  idx = plyr::count(idx)
  # In this step we populate our pmap with our randomly selected x coordinate and it's matching y coordiate from the coord object
  for(i in 1:nrow(idx)){
    pmap[P$pathwayname,coord$x[idx$x[i]], coord$y[idx$x[i]]] = idx$freq[i]
  }
  ## Print statement
  print(paste("Processed pathway",j))
}

## Save as text file:
pmap2D = lapply(rownames(pmap),function(x) cbind(as.data.frame(rep(x,nrow(pmap[x,,]))), pmap[x,,]) )
pmap2D = do.call(rbind, pmap2D)
colnames(pmap2D)[1] = "pathway"
write.table(pmap2D, file=paste0("~/Downloads/pathwayCoordinatesStack_",cellName,".txt"), sep="\t", quote=F, row.names = F)

# Finally we produce an image of the pathway map for testing:
tmp=pmap2D[pmap2D$pathway=="Cell Cycle",]
image(as.matrix(tmp[,-1]))

