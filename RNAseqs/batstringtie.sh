#!/bin/bash
for i in $(ls -d *1A)
do 
echo $i
cd $i
for file in $(ls *M.f2.NH1.sort.bam)
do
#echo $file
prefix=$(echo $i|sed 's/\(.*\)-1A/\1/g')
#echo $prefix
f1=${prefix}-1Af2.NH1.sort.bam
f2=../${prefix}-2A/${prefix}-2Af2.NH1.sort.bam
f3=../${prefix}-3A/${prefix}-3Af2.NH1.sort.bam
#f3=../{$prefix}{-3A/{$prefix}{-3A}f2.NH1.sort.bam
#echo $f1 $f2 $f3
fm=${prefix}M
fm4gtf=$fm.f2.NH1.sort.bam
#fm=${prefix}M
short=$fm.strintie
echo $short
#echo $fm
echo "#BSUB -J $fm.STR" >$short.lsf
#echo '#PBS -l nodes=1:ppn=4' >>$short.lsf
echo '#BSUB -n 4' >>$short.lsf
echo '#BSUB -M 16G' >>$short.lsf
echo '#BSUB -R span[hosts=1]' >>$short.lsf
echo '#BSUB -q normal' >>$short.lsf
echo '#BSUB -o %J.out' >>$short.lsf
echo '#BSUB -e %J.err' >>$short.lsf
#echo "samtools merge -@ 4 $fm.f2.NH1.sort.bam $f1 $f2 $f3 " >>$short.lsf
echo "stringtie -p 4 -G ../Bnafr.gtf -o $fm.f2.NH1.sort.gtf -l $fm $fm4gtf " >>$short.lsf

done
cd ..
done
