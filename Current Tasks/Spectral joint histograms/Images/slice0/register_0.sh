#!/bin/sh
#$ -j y
#$ -cwd
#$ -V

name=slice0_
outfile=$name.out
fixed=../Al/Al_0.tif
moving=../No/No_0.tif
threads=2
transform=s

antsRegistrationSyN.sh -d 2 -f $fixed -m $moving -o $name -n $threads -t $transform >> $outfile
