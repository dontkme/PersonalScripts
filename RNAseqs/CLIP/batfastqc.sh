#!/bin/bash
for file in $(ls *.fastq.gz)
do
prefix=$(echo $file | sed 's/\(.*\)\.\.fastq\.gz/\1/g')
#echo $prefix
file1=${prefix}.R1.fq.gz
#file1=${prefix}_clean_R1.fq.gz
#echo $file1
file2=${prefix}.R2.fq.gz
short=$prefix
#echo $short
#echo "#PBS -N $short.2Bnafr " >$short.lsf
echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short.fastqc" >>$short.sbatch
echo "#SBATCH --output=$short.fastp.out" >>$short.sbatch
echo "#SBATCH --error=$short.fastp.err" >>$short.sbatch
echo "#SBATCH --partition=broadwl" >>$short.sbatch
echo "#SBATCH --ntasks=2" >>$short.sbatch
echo "#SBATCH --mem=8G" >>$short.sbatch

#echo "salmon quant -i mm10.gencode.vm25.PCtrans -l A -1 ${file1} -2 ${file2} -p 2 --validateMappings -o quants/${short}_quant" >>$short.sbatch
#echo "fastp -w 2 -i $file1 -I $file2 -o $short.R1.clean.fq.gz -O $short.R2.clean.fq.gz -h $short.fastp_report.html -j $short.fastp_report.json" >>$short.sbatch
echo "module load java">>$short.sbatch
echo "fastqc -t 2 -o ../QC $file" >>$short.sbatch
echo "sbatch $short.sbatch"
done
