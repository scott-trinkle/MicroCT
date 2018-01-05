#!/bin/sh
#$ -j y
#$ -cwd
#$ -V

name=slice240
outfile=$name.out
fixed=../Al_240.tif
moving=../No_240.tif
threads=2
transform=s

antsRegistrationSyN.sh -d 2 -f $fixed -m $moving -o $name -n $threads -t $transform >> $outfile
