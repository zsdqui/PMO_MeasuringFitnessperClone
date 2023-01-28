#!/bin/sh

# Train model 
fnet train --json /raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/3D_Imaging/allenModelTraining/params150kEpochs_trainCyto.json --gpu_ids 1
