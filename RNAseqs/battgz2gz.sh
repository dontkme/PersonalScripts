#!/bin/bash
#SBATCH --partition=bigmem2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

for tgzf in $(ls *.tgz)
do
echo $tgzf
tar -zxvf $tgzf
short=$(echo "${tgzf%.*}")
echo $short
 for file in $(tar -tf $tgzf)
 do
 echo $file 
 cat $file >>${short}
 rm $file
 done
gzip $short
done
