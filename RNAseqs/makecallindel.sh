#!/bin/bash
#OUTpbsname="callallsnp.lsf"
short="callallindel"
#bamfiles=$(ls *f2.NH1.sort.bam)
ma=($(ls *M.f2.NH1.sort.bam))
echo ${#ma[*]}
#n=${#ma[*]}-1
echo >namelist.txt
for i in ${ma[*]}
do
echo $i >>namelist.txt
name+=$(echo $i | sed 's/\(.*\)\.f2.*/\1/g')
done
#echo $name
#echo ${ma[*]}
refgenome=$(ls *.fa)
usrpwd=$(pwd)
#echo $bamfiles $usrpwd $refgenome
#echo "#PBS -N callsnpall" >$OUTpbsname
#echo "#PBS -l nodes=1:ppn=1" >>$OUTpbsname
#echo "#PBS -q batch" >>$OUTpbsname
#echo "#PBS -q test" >>$OUTpbsname
#echo "#PBS -V " >>$OUTpbsname
#echo "cd $usrpwd" >>$OUTpbsname
echo "#BSUB -J callallindel" >$short.lsf
#echo '#PBS -l nodes=1:ppn=4' >>$short.lsf
echo '#BSUB -n 1' >>$short.lsf
echo '#BSUB -M 128G' >>$short.lsf
echo '#BSUB -R span[hosts=1]' >>$short.lsf
echo '#BSUB -q normal' >>$short.lsf
echo '#BSUB -o %J.out' >>$short.lsf
echo '#BSUB -e %J.err' >>$short.lsf
echo "time" >>$short.lsf
#echo "samtools mpileup -f $refgenome ${ma[*]} |java -jar /public/home/knhu/biosoft/VarScan.v2.4.4.jar mpileup2snp --vcf-sample-list namelist.txt -output-vcf 1 >$name.allsnp.vcf" >>$short.lsf
echo "samtools mpileup -f $refgenome ${ma[*]} |java -jar /public/home/knhu/biosoft/VarScan.v2.4.4.jar mpileup2indel --vcf-sample-list namelist.txt -output-vcf 1 >Allindel.vcf" >>$short.lsf
echo "Next run bsub < $short.lsf"
