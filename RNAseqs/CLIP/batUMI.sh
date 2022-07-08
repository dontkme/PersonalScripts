#!/bin/bash
#source activate /project2/xczhang/KnHu/condaenv/CLIP
for f in $(ls *trims.c.fastq.gz)
do
short=$(echo $f|sed 's/\(.*\)\.trims\.c\.fastq\.gz/\1/')
#echo $short
#echo $f
#cutadapt -f fastq --times 1 -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#cutadapt -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/fastq2collapse.pl $f - | gzip -c > $short.trim.c.fastq.gz
ncpu=4
jobname="$short.UMI"
sbatchname="$short.UMI.sbatch"
ctkpath="/home/kaininghu/xczhang/KnHu/biosoft/CTK/ctk-1.1.4"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
#echo "#SBATCH --partition=bigmem2" >>$sbatchname
echo "#SBATCH --partition=broadwl" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=10G" >>$sbatchname

#echo "perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trim.c.tag.fastq.gz" >>$sbatchname
echo "perl ${ctkpath}/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trims.c.tag.fastq.gz" >>$sbatchname
echo "sbatch $sbatchname"
done
