#!/bin/sh
#$ -j y
#$ -cwd
#$ -V

name=slice720
outfile=$name.out
fixed=../Al/Al_720.tif
moving=../No/No_720.tif
threads=2
transform=s

antsRegistrationSyN.sh -d 2 -f $fixed -m $moving -o $name -n $threads -t $transform >> $outfile
