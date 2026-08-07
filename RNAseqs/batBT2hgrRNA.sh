#!/bin/bash
#source activate /project2/xczhang/KnHu/condaenv/CLIP
#for f in  *{trim5,polyA}.fq.gz
for f in $(ls *.fq.gz | grep -E 'trim5\.fq|polyA\.fq')
do
	short=$(echo $f|sed 's/\(.*\)\.fq\.gz/\1/')
#echo $short
#echo $f
#cutadapt -f fastq --times 1 -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#cutadapt -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/fastq2collapse.pl $f - | gzip -c > $short.trim.c.fastq.gz
ncpu=4
BT2index="/sibcb1/jizhelab1/hukaining/genome/hg.rRNA_bt2_index/hg.ribosome"
jobname="$short.rmrRNA"
sbatchname="$short.rmrRNA.sbatch"
#ctkpath="/home/kaininghu/xczhang/KnHu/biosoft/CTK/ctk-1.1.4"
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
#echo "cutadapt -j 4 --times 1 -e 0.2 --discard-untrimmed  -m 18 -M 50 -u 7 -a AAAAAAAA -o $short.trim5_polyA.fq.gz $f > $short.cutadpt.log" >>$sbatchname
echo "bowtie2 -p $ncpu -x $BT2index -U $f --un-gz $short.unmapped.fq.gz -S $short.maprRNA.sam 2>$short.align_rRNA.log > $short.cutadpt.log" >>$sbatchname
echo "zcat $short.unmapped.fq.gz |awk 'NR%4==2{print length($1)}' |sort |uniq -c >$short.unmapped.fq.summary" >>$sbatchname
echo "zcat $short.unmapped.fq.gz | awk 'END{print NR/4}' > $short.unmaprRNA.readcount.txt" >>$sbatchname
echo "samtools flagstat -@ $ncpu $short.maprRNA.sam >$short.maprRNA.sam.flagstat" >>$sbatchname
#echo "umi_tools extract \

#--extract-method=string \

#--bc-pattern=NNNNNNNNNNNNNNNNNN \

#--3prime \

#-I $short.trim.fq.gz \

#-S $short.trim.umi.fq.gz \
#--log $short.umi.log " >>$sbatchname
#echo "perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trim.c.tag.fastq.gz" >>$sbatchname
#echo "cutadapt -f fastq -j 4 -u 7  -m 18 -M 50  -o $short.trim.umi.trim5.fq.gz $short.trim.umi.fq.gz > $short.trim5.log" >>$sbatchname
#echo "perl ${ctkpath}/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trims.c.tag.fastq.gz" >>$sbatchname
echo "sbatch $sbatchname"
done
