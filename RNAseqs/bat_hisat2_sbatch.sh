#!/bin/bash
reffa=$(ls *.1.ht2 |sed 's/\(.*\)\.1\.ht2/\1/')
#for i in $(ls -d *A)
#do
#cd $i
for file in $(ls *R1.clean.fq.gz)
do
prefix=$(echo $file | sed 's/\(.*\)\.R1\.clean\.fq\.gz/\1/g')
#echo $prefix
#file1=${prefix}_1.fq.gz
file1=${prefix}.R1.clean.fq.gz
#echo $file1
#file2=${prefix}_2.fq.gz
file2=${prefix}.R2.clean.fq.gz
#Ufile1=${prefix}1.unpaired.fq
#Ufile2=${prefix}2.unpaired.fq
#echo $file2
#echo "hisat2 -x $reffa -1 $file1 -2 $file2 -S $prefix.sam" >>$short.NH1.lsf
#short=$i
short=$prefix
#short=$(echo $prefix | sed 's/\(.*\)R.*/\1/')
#echo $short
#echo "#PBS -N $short.2Bnafr " >$short.lsf
#echo "#BSUB -J $short.2Bnafr" >$short.lsf
#echo '#PBS -l nodes=1:ppn=4' >>$short.lsf
#echo '#BSUB -n 4' >>$short.lsf
#echo '#BSUB -M 16G' >>$short.lsf
#echo '#BSUB -R span[hosts=1]' >>$short.lsf
#echo '#BSUB -q normal' >>$short.lsf
#echo '#BSUB -o %J.out' >>$short.lsf
#echo '#BSUB -e %J.err' >>$short.lsf
echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short.2$reffa" >>$short.sbatch
echo "#SBATCH --output=$short.2$reffa.out" >>$short.sbatch
echo "#SBATCH --error=$short.2$reffa.err" >>$short.sbatch
echo "#SBATCH --partition=broadwl" >>$short.sbatch
echo "#SBATCH --ntasks=4" >>$short.sbatch
echo "#SBATCH --mem=16G" >>$short.sbatch

#echo '#PBS -q batch' >>$short.lsf
#echo '#PBS -V ' >>$short.lsf
#echo "cd $(pwd)" >>$short.lsf
#echo "hisat2 -p 4 -x $reffa -1 $file1 -2 $file2 -U $Ufile1 -U $Ufile2 --max-intronlen 30000 -S ${short}.sam" >>$short.lsf
echo "hisat2 -p 4 -x $reffa -1 $file1 -2 $file2 --max-intronlen 299999 --no-softclip --dta -S ${short}.sam" >>$short.sbatch
echo "samtools flagstat -@ 4 ${short}.sam >$short.flagstat" >>$short.sbatch
#echo "samtools view -b -o $short.bam $short.sam" >>$short.NH1.lsf
echo "samtools sort -@ 4 -o $short.sort.bam $short.sam " >>$short.sbatch
echo "samtools index -@ 4 $short.sort.bam" >>$short.sbatch
echo "samtools view -H $short.sam >$short.NH1.sam" >>$short.sbatch
echo "grep -P \"NH:i:1\$\" $short.sam >>$short.NH1.sam" >>$short.sbatch
echo "samtools view -@ 4 -b -o $short.NH1.bam $short.NH1.sam" >>$short.sbatch
#echo "samtools view -b -o $short.bam $short.sam" >>$short.lsf
#echo "samtools view -@ 10 -f 2 -b -o Bf2.bam  B.bam
echo "samtools flagstat -@ 4 ${short}.NH1.bam >${short}.NH1.flagstat" >>$short.sbatch
echo "samtools view -@ 4 -f 2 -b -o ${short}f2.NH1.bam $short.NH1.bam">>$short.sbatch
#samtools sort -@ 10 Bf2.bam Bf2.sort
echo "samtools sort -@ 4 ${short}f2.NH1.bam -o ${short}f2.NH1.sort.bam" >>$short.sbatch
#echo "samtools sort -@ 4 ${short}.bam -o ${short}.sort.bam" >>$short.lsf
echo "samtools index -@ 4 ${short}f2.NH1.sort.bam" >>$short.sbatch
#echo "echo \"$short RNA hisat2 $reffa  done.\" " >>$short.lsf
echo "samtools flagstat -@ 4 ${short}f2.NH1.sort.bam >${short}f2.NH1.flagstat" >>$short.sbatch
echo "echo \"$short RNA hisat2 $reffa  done.\" " >>$short.sbatch

echo "sbatch $short.sbatch"
#echo $i
#done
#cd ..
done
