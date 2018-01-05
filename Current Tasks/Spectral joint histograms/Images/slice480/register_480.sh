#!/bin/sh
#$ -j y
#$ -cwd
#$ -V

name=slice480_
outfile=$name.out
fixed=../Al/Al_480.tif
moving=../No/No_480.tif
threads=2
transform=s

antsRegistrationSyN.sh -d 2 -f $fixed -m $moving -o $name -n $threads -t $transform >> $outfile
