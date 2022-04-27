#!/bin/sh
programname=$0

usage() {
  echo "usage: $programname FoF"
  echo "  FoF:           unique FoF ID associated with this acquisition"
  exit 1
  }

if [ $# -ne 4 ] ; then
  usage
fi

## Command line arguments
FoF=$1
# FoF='FoF1_211027_fluorescent.mito'

## Cellpose path, root path and target path
run_CellPose='/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code/run_CellPose_3D_v1.py'
root='/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/'
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
