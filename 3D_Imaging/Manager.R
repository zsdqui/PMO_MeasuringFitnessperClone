setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87/")
FoF="FoF00_2021Aug03"
INDIR="H03_CellposeOutput"
OUTCORRECTED="I04_PostProcessCellposeOutput"
signal1="nucleus.t_Cells_Centers.csv"
signal2="nucleus.p_Cells_Centers.csv"

## First correct segmentation output
CorrectCellposeSegmentation(FoF,signal="nucleus.t",INDIR,OUTCORRECTED,doplot=F)
CorrectCellposeSegmentation(FoF,signal="nucleus.p",INDIR,OUTCORRECTED,doplot=F)

## Next link each predicted nucleus to its closest target nucleus
setwd(paste0(OUTCORRECTED,filesep,FoF,filesep,"Cells_center_coordinates"))
OUTLINKED=paste0(getwd(),filesep,"J05_multiSignals_Linked/",FoF,filesep)
assignCompartment2Nucleus(signal2, signal1, OUTLINKED)

## Compare each predicted to its linked target nucleus
stats=compareCells(signal1, signal2,OUTLINKED)

##Plot
minmax=quantile(unlist(stats[,1:2]), c(0,1))
par(mfrow=c(2,2))
plot(stats$nucleus.t_NumPixels,stats$nucleus.p_NumPixels,pch=20,log="xy",xlim=minmax,ylim=minmax)
hist(stats$nucleus.t_IntersectingPixels,col="cyan")
hist(stats$nucleus.p_IntersectingPixels,col="cyan")

