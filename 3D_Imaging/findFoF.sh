#!/bin/sh

date=$1
filename=$2

A01='/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/NCI-N87/A01_rawData'

f=`find $A01/*$date* -type f -name README`

FoF=`grep $filename $f`
FoF="$(dirname -- $FoF)"
FoF="$(basename -- $FoF)"

echo $FoF