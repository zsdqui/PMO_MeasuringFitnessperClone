#!/bin/sh

## Command line arguments
date=$1
organelle=$2
source=$3
channel=$4
echo $date $organelle $source $channel
## Root path and step ID
root='/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87'
stepID='A01_rawData'

## Fields of view
x=`echo $organelle | sed "s/leus//g" | sed "s/Mito//g" | sed "s/plasm//g"`;
FoFs=`ls -d *$x*`
organelle=`echo "$organelle" | awk '{print tolower($0)}'`

##Iterate through fields of view
i=0
for FoF in $FoFs; do
  ## @TODO: this doesn't work if we are working with brightfield exclusively (with no fluorescent channel)
  id="FoF"$i"_"$date"_fluorescent."$organelle

  ## Target folder
  cd $FoF
  copyto=$root/$stepID/$id
  mkdir $copyto

  ## Copy files
  f=`ls *_ch$channel".tif"`
  for x in $f; do
    y=`echo $x | sed "s/$FoF/"$organelle"."$source"/g" | sed "s/_ch"$channel"//g"`
    cp $x $copyto/$y
  done

  ## Write README file <-- @TODO: split by _ and parse
  echo $FoF >> $copyto/README
  cd ..
  i=$((i+1))
done
