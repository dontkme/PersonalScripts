#!/bin/bash
for i in $(ls *.sort.bam)
do
#OUTpbsname="$i.mdp.pbs"
OUTpbsname="$i.mdp.lsf"
usrpwd=$(pwd)
#echo $bamfiles $usrpwd $refgenome
#echo "#PBS -N $OUTpbsname.2Bnafr" >$OUTpbsname
echo "#BSUB -J $OUTpbsname.2Bnafr" >$OUTpbsname
#echo "#PBS -l nodes=1:ppn=1" >>$OUTpbsname
echo "#BSUB -n 10" >>$OUTpbsname
#echo "#PBS -q batch" >>$OUTpbsname
echo "#BSUB -R span[hosts=1]" >>$OUTpbsname
echo "#BSUB -o %J.out" >>$OUTpbsname
echo "#BSUB -e %J.err" >>$OUTpbsname
echo "#BSUB -q normal" >>$OUTpbsname
#echo "#PBS -q test" >>$OUTpbsname
#echo "#PBS -V " >>$OUTpbsname
#echo "cd $usrpwd" >>$OUTpbsname
#echo "time" >>$OUTpbsname
echo "java -jar ~/biosoft/picard.jar MarkDuplicates I=$i O=$i.mdp.bam M=$i.M REMOVE_DUPLICATES=false" >>$OUTpbsname
echo "bsub <$OUTpbsname"
#echo $i
done
