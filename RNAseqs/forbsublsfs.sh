#!/bin/bash
#for i in $(ls *.pbs|sed 's/\(.*[0-9]\)\.pbs/\1/')
for i in $(ls |grep ".*[0-9]\.lsf")
do
echo $i
bsub <$i
done
