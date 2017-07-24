#!/bin/bash
pbsfile=fortrim.pbs
echo '#PBS -N TRIM ' >>$pbsfile
echo '#PBS -l nodes=1:ppn=4' >>$pbsfile
echo '#PBS -q batch' >>$pbsfile
echo '#PBS -V ' >>$pbsfile
echo "cd $(pwd)" >>$pbsfile

for file in $(ls *R1_001.fastq)
do
#prefix=$(echo $file | sed 's/\(.*\)1\.clean\.fq/\1/g')
#R1_001.fastq
prefix=$(echo $file | sed 's/\(.*\)R1_001\.fastq/\1/g')
#echo $prefix
file1=${prefix}R1_001.fastq
#echo $file1
file2=${prefix}R2_001.fastq
short=$(echo $prefix | sed 's/\(.*\)_S.*/\1/')

echo "java -jar ~/biosoft/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 4 -phred33 $file1 $file2 ${short}R1.clean.fq ${short}R1.unpaired.fq ${short}R2.clean.fq ${short}R2.unpaired.fq ILLUMINACLIP:/public/home/knhu/biosoft/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35" >>$pbsfile

done
