#!/bin/bash 
ncpu=4 
for file in $(ls -d SRR* |grep -v fastq) 
do 
echo $file 
fasterq-dump -e $ncpu -S $file 
done 