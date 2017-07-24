#!/bin/bash
#reffa=$(ls *.fa |sed 's/\(.*\)\.fa/\1/')
reffa=$(ls *.1.ht2 |sed 's/\(.*\)\.1\.ht2/\1/')
for file in $(ls *R1.clean.fq)
do
prefix=$(echo $file | sed 's/\(.*\)1\.clean\.fq/\1/g')
#echo $prefix
file1=${prefix}1.clean.fq
#echo $file1
file2=${prefix}2.clean.fq
Ufile1=${prefix}1.unpaired.fq
Ufile2=${prefix}2.unpaired.fq
#echo $file2
#echo "hisat2 -x $reffa -1 $file1 -2 $file2 -S $prefix.sam" >>$short.NH1.pbs
short=$(echo $prefix | sed 's/\(.*\)R.*/\1/')
#echo $short
echo "#PBS -N $short.2ATH " >$short.pbs
echo '#PBS -l nodes=1:ppn=4' >>$short.pbs
echo '#PBS -q batch' >>$short.pbs
echo '#PBS -V ' >>$short.pbs
echo "cd $(pwd)" >>$short.pbs
echo "hisat2 -p 4 -x $reffa -1 $file1 -2 $file2 -U $Ufile1 -U $Ufile2 --max-intronlen 10000 -S ${short}.sam" >>$short.pbs
echo "samtools flagstat ${short}.sam >$short.flagstat" >>$short.pbs
#echo "samtools view -b -o $short.bam $short.sam" >>$short.NH1.pbs
#echo "samtools view -H $short.sam >$short.NH1.sam" >>$short.NH1.pbs
#echo "grep -P \"NH:i:1\$\" $short.sam >>$short.NH1.sam" >>$short.NH1.pbs
#echo "samtools view -b -o $short.NH1.bam $short.NH1.sam" >>$short.NH1.pbs
echo "samtools view -b -o $short.bam $short.sam" >>$short.pbs
#echo "samtools view -@ 10 -f 2 -b -o Bf2.bam  B.bam
#echo "samtools view -@ 4 -f 2 -b -o ${short}f2.NH1.bam $short.NH1.bam">>$short.NH1.pbs
#samtools sort -@ 10 Bf2.bam Bf2.sort
#echo "samtools sort -@ 4 ${short}f2.NH1.bam ${short}f2.NH1.sort" >>$short.NH1.pbs
echo "samtools sort -@ 4 ${short}.bam -o ${short}.sort.bam" >>$short.pbs
echo "samtools index ${short}.sort.bam" >>$short.pbs
echo "echo \"$short RNA hisat2 $reffa  done.\" " >>$short.pbs
done
