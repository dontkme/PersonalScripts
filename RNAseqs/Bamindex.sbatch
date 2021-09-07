#!/bin/bash
#SBATCH --job-name=Bamindex
#SBATCH --output=Bamindex.out
#SBATCH --error=Bamindex.err
#SBATCH --partition=bigmem2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G

for i in $(ls *.bam)
do 
samtools index -@ 4 $i
done
