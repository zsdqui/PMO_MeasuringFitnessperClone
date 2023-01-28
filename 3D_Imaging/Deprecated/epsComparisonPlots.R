# Examine volume of nuclei
# add patchwork to requirements
library(ggplot2)
library(patchwork)

setwd("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/C06_Stats/FoF12_211110_fluorescent.nucleus")

f=list.files(pattern="stats")
X=lapply(f, read.csv)
X_=sapply(X, function(x) x$vol_nucleus.t[x$vol_nucleus.t>0])
names(X_)=f

############
## Volume ##
############
# Log transform
# Create pdf 
pdf("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/C06_Stats/FoF12_211110_fluorescent.nucleus/volumeHist_eps.pdf", width = 10, height = 10)
p=lapply(names(X_), function(x) ggplot(data.frame(X_[[x]]), aes(X_[[x]])) 
         + geom_histogram(bins = 32) + ggtitle(x)
         + scale_x_continuous(trans = "log") 
         + theme(axis.text.x = element_text(size = 5))
         + labs(x = "Nucleus Volume"))

p[[1]]+p[[2]]+p[[3]]+p[[4]]+p[[5]]+p[[6]]+p[[7]]+p[[8]]+p[[9]]+p[[10]]
p[[11]]+p[[12]]+p[[13]]+p[[14]]+p[[15]]
dev.off()

# No log transform
# Create pdf 
pdf("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/C06_Stats/FoF12_211110_fluorescent.nucleus/volumeHist_eps_noLog.pdf", width = 10, height = 10)
p=lapply(names(X_), function(x) ggplot(data.frame(X_[[x]]), aes(X_[[x]])) 
         + geom_histogram(bins = 32) + ggtitle(x)
         + theme(axis.text.x = element_text(size = 5))
         + labs(x = "Nucleus Volume"))

p[[1]]+p[[2]]+p[[3]]+p[[4]]+p[[5]]+p[[6]]+p[[7]]+p[[8]]+p[[9]]+p[[10]]
p[[11]]+p[[12]]+p[[13]]+p[[14]]+p[[15]]
dev.off()

##########
## Area ##
##########
X_=sapply(X, function(x) x$area_nucleus.t[x$area_nucleus.t>0])
names(X_)=f

# Log transform
# Create pdf 
pdf("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/C06_Stats/FoF12_211110_fluorescent.nucleus/areaHist_eps.pdf", width = 10, height = 10)
p=lapply(names(X_), function(x) ggplot(data.frame(X_[[x]]), aes(X_[[x]])) 
         + geom_histogram(bins = 32) + ggtitle(x)
         + scale_x_continuous(trans = "log") 
         + theme(axis.text.x = element_text(size = 5))
         + labs(x = "Nucleus Area"))

p[[1]]+p[[2]]+p[[3]]+p[[4]]+p[[5]]+p[[6]]+p[[7]]+p[[8]]+p[[9]]+p[[10]]
p[[11]]+p[[12]]+p[[13]]+p[[14]]+p[[15]]
dev.off()

# No log transform
# Create pdf 
pdf("/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/C06_Stats/FoF12_211110_fluorescent.nucleus/areaHist_eps_noLog.pdf", width = 10, height = 10)
p=lapply(names(X_), function(x) ggplot(data.frame(X_[[x]]), aes(X_[[x]])) 
         + geom_histogram(bins = 32) + ggtitle(x)
         + scale_x_continuous(trans = "log") 
         + theme(axis.text.x = element_text(size = 5))
         + labs(x = "Nucleus Area"))

p[[1]]+p[[2]]+p[[3]]+p[[4]]+p[[5]]+p[[6]]+p[[7]]+p[[8]]+p[[9]]+p[[10]]
p[[11]]+p[[12]]+p[[13]]+p[[14]]+p[[15]]
dev.off()

######