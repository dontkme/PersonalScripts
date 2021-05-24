#!/bin/bash
for i in $(ls *.sbatch)
do
sbatch ./${i}
done
