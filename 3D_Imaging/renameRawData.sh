#!/bin/sh
programname=$0

usage() {
	echo "usage: $programname date organelle source channel"
	echo "  date:           acquisition date"
	echo "  organelle:      origin of fluorescence signal (accepted values: mito, nucleus, cytoplasm)"
	echo "  source:         are images to be used as signal or target (accepted values: s, t for brightfield, fluorescence respectively)"
	echo "  channel:        on which channel do we find the source (e.g. 00 for fluorescence, 01 for brightfield)"
	echo "  cellLine:      cell line (default: NCI-N87)"
	exit 1
	}

if [ $# -lt 4 ] ; then
	usage
fi

## Command line arguments
date=$1
organelle=$2
source=$3
channel=$4
cellLine='NCI-N87'
if [ "$#" -eq 5 ]; then
   cellLine=$5
fi
echo $date $organelle $source $channel $cellLine

## Root path and step ID
root="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/"$cellLine
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
