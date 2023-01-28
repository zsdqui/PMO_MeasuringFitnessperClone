#!/bin/bash

FoF=$1
# FoF='FoF1_211027_fluorescent.mito'
targetorganelle="all"
if [ "$#" -ge 2 ]; then
   targetorganelle=$2
fi
cellLine='NCI-N87'
if [ "$#" -eq 3 ]; then
   cellLine=$3
fi

<<<<<<< HEAD
run_CellPose='/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/run_CellPose_3D_v1.py'
=======
run_CellPose='/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/3D_Imaging/run_CellPose_3D_v1.py'
>>>>>>> a3a43e9 (Added .sh for train/predict and cyto params)
root="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/$cellLine/"
A03=$root/A03_allenModel/$FoF
A04=$root/A04_CellposeOutput/$FoF
mkdir $A04
## temp copy to cellpose folder
cp $A03/*.tif $A04/

## Segment each 
f=`ls $A04/*p.tif  $A04/*t.tif`
conda activate cellpose
for x in $f; do
  organelle="$(basename -- $x .tif)"
  organelle=${organelle%.*}
  organelle=`echo $organelle | sed 's/mito/mitochondria/g'`
  if [[ $targetorganelle != "all" && $organelle != $targetorganelle ]]; then
	continue;
  fi
  # python3 cellPose_getCount.py -dir $x -t $organelle -GPU 1
  echo "running "$organelle" model on "$x
  python3 $run_CellPose -imgPath $x -t $organelle -GPU 1
done

## Move all coordinates one level up
f=`find $A04/All_Cells_coordinates -type f -name *_cell_*csv`
for x in $f; do 
	mv $x $A04/All_Cells_coordinates/ ; 
done
#mv $f $A04/All_Cells_coordinates/
torm=`ls -d $A04/All_Cells_coordinates/*/`
#rm -r $torm

# Writing permissions to output FoF folder
chmod -R a+w $A04
