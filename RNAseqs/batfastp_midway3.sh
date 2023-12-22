#!/bin/bash
for file in $(ls *R1.fastq)
do
prefix=$(echo $file | sed 's/\(.*\)\.R1\.fastq/\1/g')
#echo $prefix
file1=${prefix}.R1.fastq
#file1=${prefix}_clean_R1.fq.gz
#echo $file1
file2=${prefix}.R2.fastq
short=$prefix
#echo $short
#echo "#PBS -N $short.2Bnafr " >$short.lsf
echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short.fastp" >>$short.sbatch
echo "#SBATCH --output=$short.fastp.out" >>$short.sbatch
echo "#SBATCH --error=$short.fastp.err" >>$short.sbatch
echo "#SBATCH --account=pi-xczhang" >>$short.sbatch
#echo "#SBATCH --partition=broadwl" >>$short.sbatch
echo "#SBATCH --partition=caslake" >>$short.sbatch
echo "#SBATCH --ntasks=1" >>$short.sbatch
echo "#SBATCH --cpus-per-task=2" >>$short.sbatch
echo "#SBATCH --mem=8G" >>$short.sbatch

#echo "salmon quant -i mm10.gencode.vm25.PCtrans -l A -1 ${file1} -2 ${file2} -p 2 --validateMappings -o quants/${short}_quant" >>$short.sbatch
echo "fastp -w 2 -i $file1 -I $file2 -o $short.R1.clean.fq.gz -O $short.R2.clean.fq.gz -h $short.fastp_report.html -j $short.fastp_report.json" >>$short.sbatch
echo "sbatch $short.sbatch"
done
