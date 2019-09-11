#!/bin/bash

#for i in $(ls *.pbs|sed 's/\(.*[0-9]\)\.pbs/\1/')
for j in $(ls -d *A)
do
echo $j
cd $j
for i in $(ls |grep "lsf")
do
echo $i
bsub <$i
done
cd ..
done
