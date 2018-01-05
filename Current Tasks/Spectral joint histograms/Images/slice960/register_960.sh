#!/bin/sh
#$ -j y
#$ -cwd
#$ -V

name=slice960
outfile=$name.out
fixed=../Al/Al_960.tif
moving=../No/No_960.tif
threads=2
transform=s

antsRegistrationSyN.sh -d 2 -f $fixed -m $moving -o $name -n $threads -t $transform >> $outfile
