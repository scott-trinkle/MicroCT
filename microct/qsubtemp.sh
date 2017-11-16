#!/bin/sh
#$ -j y
#$ -cwd
infile=sampleall.py
runhost="$(hostname | cut -f1 -d.)"
rundate="$(date '+ %m-%d-%y %H:%M')"
outfile=runsample.out
echo "Batch job for $infile started on $runhost on $rundate" > $outfile
python < $infile >> $outfile
