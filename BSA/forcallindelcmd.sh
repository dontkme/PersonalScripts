#!/bin/bash
OUTpbsname="callallindel.pbs"
#bamfiles=$(ls *f2.NH1.sort.bam)
ma=($(ls *f2.NH1.sort.bam))
echo ${#ma[*]}
#n=${#ma[*]}-1
for i in ${ma[*]}
do
name+=$(echo $i | sed 's/\(.*\)f2.*/\1/g')
done
#echo $name
#echo ${ma[*]}
refgenome=$(ls *.fa)
usrpwd=$(pwd)
#echo $bamfiles $usrpwd $refgenome
echo "#PBS -N callindelall" >$OUTpbsname
echo "#PBS -l nodes=1:ppn=1" >>$OUTpbsname
#echo "#PBS -q batch$" >>$OUTpbsname
echo "#PBS -q test" >>$OUTpbsname
echo "#PBS -V " >>$OUTpbsname
echo "cd $usrpwd" >>$OUTpbsname
echo "time" >>$OUTpbsname
echo "samtools mpileup -f $refgenome ${ma[*]} |java -jar /public/home/knhu/biosoft/VarScan.v2.3.7.jar mpileup2indel -output-vcf 1 >$name.allindel.vcf" >>$OUTpbsname
echo "Next run qsub $OUTpbsname"
