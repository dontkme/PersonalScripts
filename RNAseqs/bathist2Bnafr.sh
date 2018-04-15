#!/bin/bash
reffa=$(ls *.1.ht2 |sed 's/\(.*\)\.1\.ht2/\1/')
for i in $(ls -d *A)
do
cd $i
for file in $(ls *_1.fq.gz)
do
prefix=$(echo $file | sed 's/\(.*\)_1\.fq\.gz/\1/g')
#echo $prefix
file1=${prefix}_1.fq.gz
#echo $file1
file2=${prefix}_2.fq.gz
#Ufile1=${prefix}1.unpaired.fq
#Ufile2=${prefix}2.unpaired.fq
#echo $file2
#echo "hisat2 -x $reffa -1 $file1 -2 $file2 -S $prefix.sam" >>$short.NH1.pbs
short=$i
#short=$(echo $prefix | sed 's/\(.*\)R.*/\1/')
#echo $short
echo "#PBS -N $short.2Bnafr " >$short.pbs
echo '#PBS -l nodes=1:ppn=4' >>$short.pbs
echo '#PBS -q batch' >>$short.pbs
echo '#PBS -V ' >>$short.pbs
echo "cd $(pwd)" >>$short.pbs
#echo "hisat2 -p 4 -x $reffa -1 $file1 -2 $file2 -U $Ufile1 -U $Ufile2 --max-intronlen 30000 -S ${short}.sam" >>$short.pbs
echo "hisat2 -p 4 -x ../$reffa -1 $file1 -2 $file2 --max-intronlen 30000 -S ${short}.sam" >>$short.pbs
echo "samtools flagstat ${short}.sam >$short.flagstat" >>$short.pbs
#echo "samtools view -b -o $short.bam $short.sam" >>$short.NH1.pbs
echo "samtools view -H $short.sam >$short.NH1.sam" >>$short.pbs
echo "grep -P \"NH:i:1\$\" $short.sam >>$short.NH1.sam" >>$short.pbs
echo "samtools view -b -o $short.NH1.bam $short.NH1.sam" >>$short.pbs
#echo "samtools view -b -o $short.bam $short.sam" >>$short.pbs
#echo "samtools view -@ 10 -f 2 -b -o Bf2.bam  B.bam
echo "samtools view -@ 4 -f 2 -b -o ${short}f2.NH1.bam $short.NH1.bam">>$short.pbs
#samtools sort -@ 10 Bf2.bam Bf2.sort
echo "samtools sort -@ 4 ${short}f2.NH1.bam -o ${short}f2.NH1.sort.bam" >>$short.pbs
#echo "samtools sort -@ 4 ${short}.bam -o ${short}.sort.bam" >>$short.pbs
echo "samtools index ${short}f2.NH1.sort.bam" >>$short.pbs
echo "echo \"$short RNA hisat2 $reffa  done.\" " >>$short.pbs

#echo $i
done
cd ..
done
