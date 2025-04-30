setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCLs/3Dbrightfield/NCI-N87/I08_3DCellProfiler_FUCCI/FoFX_2410_fluorescent.nucleus")
library(flowCore)
library(Biobase)
fucci <- read.csv(file = "object.csv")
fucci$FileName_bright=gsub(".ome.tif","",gsub("stk_0001_","",fucci$FileName_bright))
colnames(fucci)=gsub("Location_Center_","", colnames(fucci))
fucci$FileName_bright=gsub("_ch1","",fucci$FileName_bright)
fucci$FileName_bright=gsub("_ch01","",fucci$FileName_bright)
colnames(fucci) = gsub("fluor_1","green",colnames(fucci)) ## fluor_1 = green
colnames(fucci) = gsub("fluor_2","red",colnames(fucci)) ## fluor_2 = red

## Retain cell and FoF IDs in FlowJo for mapping back later
tmp=strsplit(fucci$FileName_bright,"_")
fucci$FoF =as.numeric( gsub("FoF","",sapply(tmp,"[[",1)))
fucci$Date=as.numeric(sapply(tmp,"[[",2))
fucci$CellID=1:nrow(fucci)
rownames(fucci)=as.character(fucci$CellID)

## Save for input to FlowJo
fucci_=fucci[,c("CellID","FoF","Date",grep("MeanIntensity_",colnames(fucci),value=T))]
metadata <- data.frame(
  name = colnames(fucci_),
  desc = colnames(fucci_),
  range = apply(fucci_, 2, max),
  minRange = apply(fucci_, 2, min),
  maxRange = apply(fucci_, 2, max)
)
data_ff <- flowFrame(as.matrix(fucci_), parameters = AnnotatedDataFrame(metadata))
write.FCS(data_ff, "~/Downloads/2410_fluorescent.nucleus_fucci.fcs")
plot(fucci_$Intensity_MeanIntensity_green,fucci_$Intensity_MeanIntensity_red, log='xy', pch=20)

## plot and save FlowJo output 'concat_1.csv':
cc=read.csv("~/Desktop/concat_1.csv")
fucci=cbind(fucci[as.character(cc$CellID),],as.matrix(cc$SampleID))
colnames(fucci)[ncol(fucci)] = "cellCycle"
fucci$cellCycle=gsub("1","S",fucci$cellCycle)
fucci$cellCycle=gsub("4","G1",fucci$cellCycle)
fucci$cellCycle=gsub("2","G1S",fucci$cellCycle)
fucci$cellCycle=gsub("3","G2M",fucci$cellCycle)

col=1:4
names(col)=unique(fucci$cellCycle)
plot(fucci$Intensity_MeanIntensity_green,fucci$Intensity_MeanIntensity_red, log='xy', pch=20, col=col[fucci$cellCycle])
legend("topright",names(col),fill=col)

fr=plyr::count(fucci$cellCycle)
fr$freq=fr$freq/sum(fr$freq)
write.csv(fucci,file = "object_cellCycle.csv",row.names = F, quote = F)
