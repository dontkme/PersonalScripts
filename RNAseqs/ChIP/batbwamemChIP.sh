#!/bin/bash
picardpath="/project2/xczhang/KnHu/bin/picard.jar"
refbwa="/project2/xczhang/genomes/GENCODE_M25_BWA_index/GRCm38.p6.genome.fa"
ncpu=4
for i in $(ls *.clean.fq.gz)
do
short=$(echo $i|sed 's/\(.*\)\.clean.fq.gz/\1/')
#echo $short
jobname="$short.BWA.aln"
sbatchname="$short.BWA.aln.sbatch"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
#echo "#SBATCH --partition=bigmem2" >>$sbatchname
echo "#SBATCH --partition=broadwl" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=20G" >>$sbatchname

f=$short.aln
#echo "bwa mem -M -t $ncpu $refbwa $i | gzip -c > $f.sam.gz" >>$sbatchname
echo "bwa aln -t $ncpu -n 0.06 -q 20 $refbwa $i > $f.sai">>$sbatchname

#echo "bwa samse /project2/xczhang/genomes/GENCODE_M25_BWA_index/GRCm38.p6.genome.fa $f.sai $f.trims.c.tag.fastq.gz | gzip -c > $f.sam.gz" >>$sbatchname
echo "bwa samse $refbwa $f.sai $i | gzip -c > $f.sam.gz" >>$sbatchname

echo "samtools view -b -@ $ncpu -o $f.bam $f.sam.gz" >>$sbatchname

echo "samtools sort -o $f.sort.bam $f.bam" >>$sbatchname
echo "samtools index $f.sort.bam">>$sbatchname
echo "samtools flagstat $f.sort.bam >$f.flagstat">>$sbatchname

echo "samtools view -@ $ncpu -bF 4 $f.sort.bam > $f.F4.sort.bam" >>$sbatchname

echo "samtools index $f.F4.sort.bam">>$sbatchname

echo "module load java" >>$sbatchname
echo "java -jar $picardpath MarkDuplicates I=$f.F4.sort.bam O=$f.rmdup.bam M=$f.rmdup.metrics.txt REMOVE_DUPLICATES=true">>$sbatchname
echo "samtools sort -@ $ncpu -o $f.rmdup.sort.bam $f.rmdup.bam" >>$sbatchname
echo "samtools index $f.rmdup.sort.bam">>$sbatchname 
echo "samtools flagstat $f.rmdup.sort.bam >$f.rmdup.flagstat">>$sbatchname
echo "sbatch $sbatchname"
done
