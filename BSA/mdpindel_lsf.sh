#!/bin/bash
OUTpbsname="callallindel.lsf"
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
#echo "#PBS -N callindelall" >$OUTpbsname
#echo "#PBS -l nodes=1:ppn=1" >>$OUTpbsname
#echo "#PBS -q batch" >>$OUTpbsname
#echo "#PBS -q test" >>$OUTpbsname
#echo "#PBS -V " >>$OUTpbsname
echo "#BSUB -J callindelall" >$OUTpbsname
echo "#BSUB -n 10" >>$OUTpbsname
echo "#BSUB -R span[hosts=1]" >>$OUTpbsname
echo "#BSUB -o %J.out" >>$OUTpbsname
echo "#BSUB -e %J.err" >>$OUTpbsname
echo "#BSUB -q normal" >>$OUTpbsname
#echo "cd $usrpwd" >>$OUTpbsname
#echo "time" >>$OUTpbsname
echo "samtools mpileup -f $refgenome ${ma[*]} |java -jar /public/home/knhu/biosoft/VarScan.v2.4.4.jar mpileup2indel --vcf-sample-list namelist.txt --output-vcf 1 >$name.allindel.vcf" >>$OUTpbsname
echo "bsub <$OUTpbsname"
