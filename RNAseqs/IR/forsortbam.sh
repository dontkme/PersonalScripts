#!/bin/bash
ncpu=4
for i in $(ls -d *rep*.od)
do
cd $i
samtools sort -@ $ncpu -o $i.sort.bam Unsorted.bam
samtools index $i.sort.bam
cd ..
done
