#for i in $(ls -d *A)
#do
#cd $i
for file in $(ls *S*_R1_001.fastq.gz)
do
prefix=$(echo $file | sed 's/\(.*\)_R1_001\.fastq\.gz/\1/g')
echo $prefix
file1=${prefix}_R1_001.fastq.gz
#file1=${prefix}_clean_R1.fq.gz
#echo $file1
file2=${prefix}_R2_001.fastq.gz
#file2=${prefix}_clean_R2.fq.gz
#Ufile1=${prefix}1.unpaired.fq
#Ufile2=${prefix}2.unpaired.fq
#echo $file2
#echo "hisat2 -x $reffa -1 $file1 -2 $file2 -S $prefix.sam" >>$short.NH1.lsf
#short=$i
#short=$prefix
short=$(echo $prefix | sed 's/.*\(S.*\)/\1/')
echo $short
#echo "#PBS -N $short.2Bnafr " >$short.lsf
echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short.2mm10" >>$short.sbatch
echo "#SBATCH --output=$short.out" >>$short.sbatch
echo "#SBATCH --error=$short.err" >>$short.sbatch
echo "#SBATCH --partition=broadwl" >>$short.sbatch
echo "#SBATCH --ntasks=2" >>$short.sbatch
echo "#SBATCH --mem=8G" >>$short.sbatch

#echo "source activate salmon" >>$short.sbatch
echo "module load java" >>$short.sbatch
#echo "salmon quant -i mm10.gencode.vm25.PCtrans -l A -1 ${file1} -2 ${file2} -p 2 --validateMappings -o quants/${short}_quant" >>$short.sbatch
echo "java -jar -Xmx8G ~/xczhang/KnHu/biosoft/picard.jar FastqToSam -F1 $file1 -F2 $file2 -O $short.unmapped.bam -SM $short -RG $short" >>$short.sbatch


#echo "#BSUB -J $short.2Bnafr" >$short.lsf
#echo '#PBS -l nodes=1:ppn=4' >>$short.lsf
#echo '#BSUB -n 4' >>$short.lsf
#echo '#BSUB -M 16G' >>$short.lsf
#echo '#BSUB -R span[hosts=1]' >>$short.lsf
#echo '#BSUB -q normal' >>$short.lsf
#echo '#BSUB -o %J.out' >>$short.lsf
#echo '#BSUB -e %J.err' >>$short.lsf
#echo '#PBS -q batch' >>$short.lsf
#echo '#PBS -V ' >>$short.lsf
#echo "cd $(pwd)" >>$short.lsf
#echo "hisat2 -p 4 -x $reffa -1 $file1 -2 $file2 -U $Ufile1 -U $Ufile2 --max-intronlen 30000 -S ${short}.sam" >>$short.lsf

#echo "hisat2 -p 4 -x ../$reffa -1 $file1 -2 $file2 --max-intronlen 20000 --dta -S ${short}.sam" >>$short.lsf
#echo "samtools flagstat -@ 4 ${short}.sam >$short.flagstat" >>$short.lsf
#echo "samtools view -b -o $short.bam $short.sam" >>$short.NH1.lsf
#echo "samtools sort -@ 4 -o $short.sort.bam $short.sam " >>$short.lsf
#echo "samtools index -@ 4 $short.sort.bam" >>$short.lsf
#echo "samtools view -H $short.sam >$short.NH1.sam" >>$short.lsf
#echo "grep -P \"NH:i:1\$\" $short.sam >>$short.NH1.sam" >>$short.lsf
#echo "samtools view -@ 4 -b -o $short.NH1.bam $short.NH1.sam" >>$short.lsf
#echo "samtools view -b -o $short.bam $short.sam" >>$short.lsf
#echo "samtools view -@ 10 -f 2 -b -o Bf2.bam  B.bam
#echo "samtools flagstat -@ 4 ${short}.NH1.bam >${short}.NH1.flagstat" >>$short.lsf
#echo "samtools view -@ 4 -f 2 -b -o ${short}f2.NH1.bam $short.NH1.bam">>$short.lsf
#samtools sort -@ 10 Bf2.bam Bf2.sort
#echo "samtools sort -@ 4 ${short}f2.NH1.bam -o ${short}f2.NH1.sort.bam" >>$short.lsf
#echo "samtools sort -@ 4 ${short}.bam -o ${short}.sort.bam" >>$short.lsf
#echo "samtools index -@ 4 ${short}f2.NH1.sort.bam" >>$short.lsf
#echo "echo \"$short RNA hisat2 $reffa  done.\" " >>$short.lsf
#echo "samtools flagstat -@ 4 ${short}f2.NH1.sort.bam >${short}f2.NH1.flagstat" >>$short.lsf
#echo "echo \"$short RNA hisat2 $reffa  done.\" " >>$short.lsf

#echo $i
done
#cd ..
#done
