#!/bin/bash
#$ -t 1-13 -m a -cwd -N CLIP
for f in Fox1_1 Fox2_1 Fox3_1
do
#f=${files[$SGE_TASK_ID-1]}
ncpu=4
jobname="$f.bwa"
sbatchname="$f.bwa.sbatch"
ctkpath="/home/kaininghu/xczhang/KnHu/biosoft/CTK/ctk-1.1.4"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
echo "#SBATCH --partition=bigmem2" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=60G" >>$sbatchname

echo "bwa aln -t $ncpu -n 0.06 -q 20 /project2/xczhang/genomes/GENCODE_M25_BWA_index/GRCm38.p6.genome.fa $f.trim.c.tag.fastq.gz > $f.sai">>$sbatchname

echo "bwa samse /project2/xczhang/genomes/GENCODE_M25_BWA_index/GRCm38.p6.genome.fa $f.sai $f.trim.c.tag.fastq.gz | gzip -c > $f.sam.gz" >>$sbatchname

echo "perl /home/kaininghu/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file $f.mutation.txt $f.sam.gz $f.tag.bed">>$sbatchname

echo "perl ${ctkpath}/tag2collapse.pl -v -big --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name $f.tag.bed $f.tag.uniq.bed">>$sbatchname

echo "awk '{print \$3-\$2}' $f.tag.uniq.bed | sort -n | uniq -c | awk '{print \$2\"\t\"\$1}' > $f.uniq.len.dist.txt">>$sbatchname

echo "python ${ctkpath}/joinWrapper.py $f.mutation.txt $f.tag.uniq.bed 4 4 N $f.tag.uniq.mutation.txt">>$sbatchname

echo "perl ${ctkpath}/bed2rgb.pl -v -col "128,0,0" $f.tag.uniq.bed $f.tag.uniq.rgb.bed">>$sbatchname
#### need modify -col


echo "sbatch $sbatchname"
done
