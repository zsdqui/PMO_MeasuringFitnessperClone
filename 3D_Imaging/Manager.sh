code='/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/code'

# ## Renome data
# cd /raid/crdlab/ix2/andrew/PMO_MeasuringFitnessperClone/3D_Imaging/Project_Nucleus/Data_Images/210803_Nuclear/Images/
# renameRawData=$code/renameRawData.sh
# date='210803'
# organelle='Nucleus'
# source='t'
# channel='00' #fluorescense
# $renameRawData $date $organelle $source $channel
# source='s'
# channel='01' #brightfield
# $renameRawData $date $organelle $source $channel

##################
####Single FoF####
# FoF='FoF5_210803_fluorescent.nucleus'
FoF='FoF11_211215_fluorescent.mito'
# FoF='FoF13_220228_fluorescent.cytoplasm'

## @TODO ome.tiff conversion

## AllenModel
$code/allenModelPredict.sh $FoF nucleus &&
$code/allenModelPredict.sh $FoF mito 

## CellPose
$code/cellPoseSegment.sh $FoF

##################
#### All FoFs ####
cd /raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/A02_ometiffconversion

## ome.tiff conversion
# python $code/stack2ometiff_conversion.py
f=`ls -d *cyto*ome.tif`

## AllenModel
for x in $f; do 
  FoF="$(basename -- $x .ome.tif)"
  $code/allenModelPredict.sh $FoF nucleus
  $code/cellPoseSegment.sh $FoF
  # $code/allenModelPredict.sh $FoF mito
done

# ## CellPose
# for x in $f; do 
#   FoF="$(basename -- $x .ome.tif)"
#   $code/cellPoseSegment.sh $FoF
# done
