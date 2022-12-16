#!/bin/sh
programname=$0

usage() {
	echo "usage: $programname date filename"
	echo "  date:           acquisition date"
	echo "  filename:      unique FoF ID associated with this acquisition"
	echo "  cellLine:      cell line (default: NCI-N87)"
	exit 1
	}

if [ $# -lt 2 ] ; then
	usage
fi

## Command line arguments
date=$1
filename=$2
cellLine='NCI-N87'
if [ "$#" -eq 3 ]; then
   cellLine=$3
fi

## Root path 
A01="/raid/crdlab/ix1/Projects/M005_MeasuringFitnessPerClone_2019/data/GastricCancerCLs/3Dbrightfield/$cellLine/A01_rawData"

# Find the filename by FoF and match with directory and basename
f=`find $A01/*$date* -type f -name README`

FoF=`grep $filename'$' $f`
FoF="$(dirname -- $FoF)"
FoF="$(basename -- $FoF)"

echo $FoF
