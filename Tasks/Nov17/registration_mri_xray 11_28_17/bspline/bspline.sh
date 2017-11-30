#!/bin/sh
#$ -j y
#$ -cwd
runhost="$(hostname | cut -f1 -d.)"
rundate="$(date '+ %m/%d/%y %H:%M')"

regtype=bspline
outfile=$regtype.out
mri=../XRAY_MRI_DATA/MRI_structural.nii
xray=../XRAY_MRI_DATA/xray_small.nii
threads=8
transform=b

echo "Batch job for $infile started on $runhost on $rundate" > $outfile
/home/trinkle/bin/ants/bin/antsRegistrationSyN.sh -d 3 -f $mri -m $xray -o $regtype -n $threads -t $transform >> $outfile
