#!/bin/sh

FoF=$1
organelle=$2

A02="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/A02_ometiffconversion"
A03="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/A03_allenModel/"$FoF
json=$A02/$FoF.$organelle".json"
csv=$A02/$FoF".csv"
target=${FoF##*.}

# Prepare json and csv file
cp /raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/allenModel/Template.json $json
cp /raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/allenModel/Template.csv $csv
sed -i "s/organelleXX/$organelle/g" $json
sed -i "s/FoFXX/$FoF/g" $json
sed -i "s/FoFXX/$FoF/g" $csv

# Run the predict function. The test images path must be present in $csv file. 
fnet predict --json $json --gpu_ids 1

# Rename Allen Model output
mv $A03/tifs/0_prediction_c0.$organelle".p.tif" $A03/$organelle".p.tif" 
# Only if input includes fluorescence signal:
if [ "$FoF" != "$target" ]; then
	mv $A03/tifs/0_target.tif $A03/$target".t.tif"
fi

# Writing permissions to output FoF folder
chmod -R a+w $A03
