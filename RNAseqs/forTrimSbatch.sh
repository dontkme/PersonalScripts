#!/bin/bash

for file in $(ls *.fastq.gz)
do
#prefix=$(echo $file | sed 's/\(.*\)1\.clean\.fq/\1/g')
#R1_001.fastq
prefix=$(echo $file | sed 's/\(.*\)\.fastq\.gz/\1/g')
#echo $prefix
pbsfile=${prefix}.trim.sbatch
file1=${prefix}.fastq.gz
cleanfile=${prefix}.clean.fq.gz
#echo $file1
file2=${prefix}R2_001.fastq
#short=$(echo $prefix | sed 's/\(.*\)_S.*/\1/')
short=${prefix}.trim
#echo $short

echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short" >>$short.sbatch
echo "#SBATCH --output=$short.out" >>$short.sbatch
echo "#SBATCH --error=$short.err" >>$short.sbatch
echo "#SBATCH --partition=broadwl" >>$short.sbatch
echo "#SBATCH --ntasks=2" >>$short.sbatch
echo "#SBATCH --mem=16G" >>$short.sbatch

echo "module load java">>$short.sbatch
echo "java -jar /home/kaininghu/xczhang/KnHu/biosoft/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 2 -phred33 $file1 $cleanfile ILLUMINACLIP:/home/kaininghu/xczhang/KnHu/biosoft/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" >>$short.sbatch
echo "sbatch $short.sbatch"

done
