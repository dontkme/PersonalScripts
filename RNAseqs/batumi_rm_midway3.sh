#!/bin/bash
for file in $(ls *_L007_R1_001.fastq.gz)
do
prefix=$(echo $file | sed 's/\(.*\)_L007_R1.*/\1/g')
#echo $prefix
file1=${prefix}_L007_R1_001.fastq.gz
#file1=${prefix}_clean_R1.fq.gz
#echo $file1
file2=${prefix}.R2.fastq
short=$prefix
#echo $short
#echo "#PBS -N $short.2Bnafr " >$short.lsf
echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short.rmumi" >>$short.sbatch
echo "#SBATCH --output=$short.rmumi.out" >>$short.sbatch
echo "#SBATCH --error=$short.rmumi.err" >>$short.sbatch
echo "#SBATCH --account=pi-xczhang" >>$short.sbatch
#echo "#SBATCH --partition=broadwl" >>$short.sbatch
echo "#SBATCH --partition=caslake" >>$short.sbatch
echo "#SBATCH --ntasks=1" >>$short.sbatch
echo "#SBATCH --cpus-per-task=2" >>$short.sbatch
echo "#SBATCH --mem=16G" >>$short.sbatch

#echo "salmon quant -i mm10.gencode.vm25.PCtrans -l A -1 ${file1} -2 ${file2} -p 2 --validateMappings -o quants/${short}_quant" >>$short.sbatch
#echo "fastp -w 2 -i $file1 -I $file2 -o $short.R1.clean.fq.gz -O $short.R2.clean.fq.gz -h $short.fastp_report.html -j $short.fastp_report.json" >>$short.sbatch
#echo "fastp -w 2 -i $file1 -o $short.R1.clean.fq.gz -h $short.fastp_report.html -j $short.fastp_report.json" >>$short.sbatch
echo "umi_tools extract --extract-method=regex -p '.+(?P<umi_1>.{8})(?P<discard_1>AGATCGGAAGAGCACA){s<=2}(?P<discard_2>.*)$' -I $file1 -S $short.rmumi.fastq.gz " >>$short.sbatch
echo "sbatch $short.sbatch"
done
