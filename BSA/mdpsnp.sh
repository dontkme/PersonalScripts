#!/bin/bash
OUTpbsname="callallsnp.pbs"
#bamfiles=$(ls *f2.NH1.sort.bam)
ma=($(ls *f2.NH1.sort.bam.mdp.bam))
echo ${#ma[*]}
#n=${#ma[*]}-1
echo >namelist.txt
for i in ${ma[*]}
do
echo $i >>namelist.txt
name+=$(echo $i | sed 's/\(.*\)f2.*/\1/g')
done
#echo $name
#echo ${ma[*]}
refgenome=$(ls *.fa)
usrpwd=$(pwd)
#echo $bamfiles $usrpwd $refgenome
echo "#PBS -N callsnpall" >$OUTpbsname
echo "#PBS -l nodes=1:ppn=1" >>$OUTpbsname
#echo "#PBS -q batch" >>$OUTpbsname
echo "#PBS -q test" >>$OUTpbsname
echo "#PBS -V " >>$OUTpbsname
echo "cd $usrpwd" >>$OUTpbsname
echo "time" >>$OUTpbsname
echo "samtools mpileup -f $refgenome ${ma[*]} |java -jar /public/home/knhu/biosoft/VarScan.v2.4.2.jar mpileup2snp --vcf-sample-list namelist.txt -output-vcf 1 >$name.allsnp.vcf" >>$OUTpbsname
echo "Next run qsub $OUTpbsname"
