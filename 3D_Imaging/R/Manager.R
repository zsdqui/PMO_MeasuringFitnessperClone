# setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87")
# conda activate r_env
setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/R")
source("CorrectCellposeSegmentation.R")
source("assignCompartment2Nucleus.R")
source("compareCells.R")
setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87")
library(matlab)

## Input and output:
FoF="FoF0_211007_fluorescent.nucleus"
INDIR="A04_CellposeOutput"
OUTCORRECTED="A05_PostProcessCellposeOutput"
OUTLINKED="A06_multiSignals_Linked"
signal1="nucleus.t_Cells_Centers.csv"
signal2="nucleus.p_Cells_Centers.csv"

## First correct segmentation output
dir.create(OUTCORRECTED)
CorrectCellposeSegmentation(FoF,signal="nucleus.t",INDIR,OUTCORRECTED,doplot=F)
CorrectCellposeSegmentation(FoF,signal="nucleus.p",INDIR,OUTCORRECTED,doplot=F)

## Next link each predicted nucleus to its closest target nucleus
OUTLINKED=paste0(getwd(),filesep,OUTLINKED,filesep,FoF,filesep)
setwd(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates"))
assignCompartment2Nucleus(signal2, signal1, OUTLINKED)

## Compare each predicted to its linked target nucleus
stats=compareCells(signal1, signal2,OUTLINKED)

##Plot
minmax=quantile(unlist(stats[,1:2]), c(0,1))
par(mfrow=c(2,2))
plot(stats$nucleus.t_NumPixels,stats$nucleus.p_NumPixels,pch=20,log="xy",xlim=minmax,ylim=minmax)
hist(stats$nucleus.t_IntersectingPixels,col="cyan")
hist(stats$nucleus.p_IntersectingPixels,col="cyan")

