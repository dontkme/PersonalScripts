#!/bin/bash
#for i in $(ls *.pbs|sed 's/\(.*[0-9]\)\.pbs/\1/')
for j in $(ls -d *A)
do 
cd $j
echo $j
for i in $(ls |grep ".*\.pbs")
do
echo $i
qsub $i
done
cd ..
done
