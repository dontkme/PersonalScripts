#!/bin/bash
for i in $(ls *.sort.bam)
do
OUTpbsname="$i.mdp.pbs"
usrpwd=$(pwd)
#echo $bamfiles $usrpwd $refgenome
echo "#PBS -N $OUTpbsname.2Bnafr" >$OUTpbsname
echo "#PBS -l nodes=1:ppn=1" >>$OUTpbsname
echo "#PBS -q batch" >>$OUTpbsname
#echo "#PBS -q test" >>$OUTpbsname
echo "#PBS -V " >>$OUTpbsname
echo "cd $usrpwd" >>$OUTpbsname
echo "time" >>$OUTpbsname
echo "java -jar ~/biosoft/picard.jar MarkDuplicates I=$i O=$i.mdp.bam M=$i.M REMOVE_DUPLICATES=false" >>$OUTpbsname
echo "Next run qsub $OUTpbsname"
#echo $i
done
