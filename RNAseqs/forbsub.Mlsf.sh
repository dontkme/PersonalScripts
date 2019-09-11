#!/bin/bash

#for i in $(ls *.pbs|sed 's/\(.*[0-9]\)\.pbs/\1/')
for j in $(ls -d *1A)
do
echo $j
cd $j
for i in $(ls |grep "merge.lsf")
do
echo $i
bsub <$i
done
cd ..
done
