#!/bin/bash
#source activate /project2/xczhang/KnHu/condaenv/CLIP
for f in $(ls *-Rp_1.fq.gz)
do
	short=$(echo $f|sed 's/\(.*-Rp_1\)\.fq\.gz/\1/')
#echo $short
#echo $f
#cutadapt -f fastq --times 1 -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#cutadapt -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/fastq2collapse.pl $f - | gzip -c > $short.trim.c.fastq.gz
ncpu=4
jobname="$short.cutAdptUMI"
sbatchname="$short.cutAdptUMI.sbatch"
ctkpath="/home/kaininghu/xczhang/KnHu/biosoft/CTK/ctk-1.1.4"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
#echo "#SBATCH --partition=bigmem2" >>$sbatchname
echo "#SBATCH --partition=default1" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=20G" >>$sbatchname

echo "source ~/.bashrc">>$sbatchname
echo "mamba activate riborf">>$sbatchname
echo "cutadapt  -j 4 --times 1 -O 5 -e 0.2 --discard-untrimmed  -m 32 -M 75 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -o $short.trim.fq.gz $f > $short.cutadpt.log" >>$sbatchname
echo "zcat $short.trim.fq.gz |awk 'NR%4==2{print length($1)}' |sort |uniq -c >$short.trim.fq.summary" >>$sbatchname

echo "umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNNNNNNNNNN --3prime -I $short.trim.fq.gz -S $short.trim.umi.fq.gz --log $short.umi.log " >>$sbatchname
#echo "perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trim.c.tag.fastq.gz" >>$sbatchname
echo "zcat $short.trim.umi.fq.gz |awk 'NR%4==2{print length($1)}' |sort |uniq -c >$short.trim.umi.fq.summary" >>$sbatchname
echo "cutadapt -j 4 -u 7  -m 18 -M 50  -o $short.trim.umi.trim5.fq.gz $short.trim.umi.fq.gz > $short.trim5.log" >>$sbatchname
#echo "perl ${ctkpath}/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trims.c.tag.fastq.gz" >>$sbatchname
echo "zcat $short.trim.umi.trim5.fq.gz |awk 'NR%4==2{print length($1)}' |sort |uniq -c >$short.trim.umi.trim5.fq.summary" >>$sbatchname
echo "sbatch $sbatchname"
done
