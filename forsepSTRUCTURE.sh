#!/bin/bash  
#echo '#PBS -N hisRNANH1 ' >>$short.NH1.pbs
#echo '#PBS -l nodes=1:ppn=4' >>$short.NH1.pbs
#echo '#PBS -q batch' >>$short.NH1.pbs
#echo '#PBS -V ' >>$short.NH1.pbs
#echo "cd $(pwd)" >>$short.NH1.pbs
#reffa=''
#reffa=$(ls *.fa |sed 's/\(.*\)\.fa/\1/')
userID="CLL"
inputfile="CLL419-t.txt"
mainparam="CLLmainpararms"
#echo $reffa
#for file in $(ls *1.clean.fq)  
for i in {1..5}
do  
mkdir "$userID.$i"
	for j in {1..10}
	do 
	prefix=$userID.$i.K$j
#prefix=$(echo $file | sed 's/\(.*\)1\.clean\.fq/\1/g')
#echo $prefix
#file1=${prefix}1.clean.fq 
#echo $file1
#file2=${prefix}2.clean.fq
#echo $file2
#echo "hisat2 -x $reffa -1 $file1 -2 $file2 -S $prefix.sam" >>$short.NH1.pbs
short=$(echo $prefix | sed 's/\(.*\)_H.*/\1/')
#echo $short
	echo "#PBS -N $prefix " >$prefix.pbs
	echo '#PBS -l nodes=1:ppn=1' >>$prefix.pbs
	echo '#PBS -q batch' >>$prefix.pbs
	echo '#PBS -V ' >>$prefix.pbs
	echo "cd /public/home/knhu/biosoft/console" >>$prefix.pbs
	echo "./structure -K $j -o $userID.$i/$j -m $mainparam -i $inputfile" >>$prefix.pbs
#echo "cd $(pwd)" >>$short.NH1.pbs
#echo "hisat2 -p 4 -x $reffa -1 $file1 -2 $file2 -S ${short}.sam" >>$short.NH1.pbs
#echo "samtools flagstat ${short}.sam >$short.flagstat" >>$short.NH1.pbs
#echo "samtools view -b -o $short.bam $short.sam" >>$short.NH1.pbs
#echo "samtools view -H $short.sam >$short.NH1.sam" >>$short.NH1.pbs
#echo "grep -P \"NH:i:1\$\" $short.sam >>$short.NH1.sam" >>$short.NH1.pbs
#echo "samtools view -b -o $short.NH1.bam $short.NH1.sam" >>$short.NH1.pbs
#echo "samtools view -@ 10 -f 2 -b -o Bf2.bam  B.bam
#echo "samtools view -@ 4 -f 2 -b -o ${short}f2.NH1.bam $short.NH1.bam">>$short.NH1.pbs
#samtools sort -@ 10 Bf2.bam Bf2.sort
#echo "samtools sort -@ 4 ${short}f2.NH1.bam ${short}f2.NH1.sort" >>$short.NH1.pbs
#samtools index Bf2.sort.bam
#echo "samtools index ${short}f2.NH1.sort.bam" >>$short.NH1.pbs
#echo "echo \"$short RNA hisat2 $reffa NH1 done.\" " >>$short.NH1.pbs
	
	echo "qsub $prefix.pbs">>qsuball.sh
	done
done  
