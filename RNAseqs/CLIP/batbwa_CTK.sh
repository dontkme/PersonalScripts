#!/bin/bash
#source activate /project2/xczhang/KnHu/condaenv/CLIP
for file in $(ls *trims.c.tag.fastq.gz)
do
short=$(echo $file|sed 's/\(.*\)\.trims\.c\.tag\.fastq\.gz/\1/')
#echo $short
#echo $f
#cutadapt -f fastq --times 1 -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#cutadapt -e 0.1 -O 1 --quality-cutoff 5 -m 20 -a TCGTATGCCGTCTTCTGCTTG -o $short.trim.fastq.gz $f > $short.cutadpt.log
#perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/fastq2collapse.pl $f - | gzip -c > $short.trim.c.fastq.gz
ncpu=4
jobname="$short.BWA"
sbatchname="$short.BWA.sbatch"
ctkpath="/home/kaininghu/xczhang/KnHu/biosoft/CTK/ctk-1.1.4"
refbwa="/project2/xczhang/genomes/GENCODE_M25_BWA_index/GRCm38.p6.genome.fa"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
echo "#SBATCH --partition=bigmem2" >>$sbatchname
#echo "#SBATCH --partition=broadwl" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=10G" >>$sbatchname

f=$short
#echo "perl /project2/xczhang/KnHu/biosoft/CTK/ctk-1.1.4/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trim.c.tag.fastq.gz" >>$sbatchname
#echo "perl ${ctkpath}/stripBarcode.pl -format fastq -len 5  $f - | gzip -c > $short.trim.c.tag.fastq.gz" >>$sbatchname

#echo "bwa aln -t $ncpu -n 0.06 -q 20 /project2/xczhang/genomes/GENCODE_M25_BWA_index/GRCm38.p6.genome.fa $f.trims.c.tag.fastq.gz > $f.sai">>$sbatchname
echo "bwa aln -t $ncpu -n 0.06 -q 20 $refbwa $f.trims.c.tag.fastq.gz > $f.sai">>$sbatchname

#echo "bwa samse /project2/xczhang/genomes/GENCODE_M25_BWA_index/GRCm38.p6.genome.fa $f.sai $f.trims.c.tag.fastq.gz | gzip -c > $f.sam.gz" >>$sbatchname
echo "bwa samse $refbwa $f.sai $f.trims.c.tag.fastq.gz | gzip -c > $f.sam.gz" >>$sbatchname

echo "samtools view -b -@ $ncpu -o $f.bam $f.sam.gz" >>$sbatchname

echo "samtools sort -o $f.sort.bam $f.bam" >>$sbatchname
echo "samtools index $f.sort.bam">>$sbatchname
echo "samtools flagstat $f.sort.bam >$f.flagstat">>$sbatchname

echo "perl ${ctkpath}/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file $f.mutation.txt $f.sam.gz $f.tag.bed">>$sbatchname

echo "perl ${ctkpath}/tag2collapse.pl -v -big --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name $f.tag.bed $f.tag.uniq.bed">>$sbatchname

echo "awk '{print \$3-\$2}' $f.tag.uniq.bed | sort -n | uniq -c | awk '{print \$2\"\t\"\$1}' > $f.uniq.len.dist.txt">>$sbatchname

echo "python ${ctkpath}/joinWrapper.py $f.mutation.txt $f.tag.uniq.bed 4 4 N $f.tag.uniq.mutation.txt">>$sbatchname

echo "perl ${ctkpath}/bed2rgb.pl -v -col "128,0,0" $f.tag.uniq.bed $f.tag.uniq.rgb.bed">>$sbatchname
#### need modify -col
#echo "perl ${ctkpath}/bed2annotation.pl -dbkey mm10 -ss -big -region -v -summary $f.tag.uniq.annot.summary.txt $f.tag.uniq.rgb.bed $f.tag.uniq.annot.txt">>$sbatchname
echo "perl ${ctkpath}/bed2annotation.pl -dbkey hg38 -ss -big -region -v -summary $f.tag.uniq.annot.summary.txt $f.tag.uniq.rgb.bed $f.tag.uniq.annot.txt">>$sbatchname

echo "perl ${ctkpath}/tag2profile.pl -v -big -ss -exact -of bedgraph -n ″${f}_Unique_Tag_Profile″ $f.tag.uniq.rgb.bed $f.tag.uniq.bedgraph">>$sbatchname

echo "perl ${ctkpath}/tag2peak.pl -big -ss -v --valley-seeking --valley-depth 0.9 $f.tag.uniq.rgb.bed $f.tag.uniq.peak.bed --out-boundary $f.tag.uniq.peak.boundary.bed --out-half-PH $f.tag.uniq.peak.halfPH.bed">>$sbatchname

echo "perl ${ctkpath}/bed2annotation.pl -dbkey mm10 -ss -big -region -v -summary $f.tag.uniq.peak.annot.summary.txt $f.tag.uniq.peak.bed $f.tag.uniq.peak.annot.txt">>$sbatchname

echo "perl ${ctkpath}/tag2peak.pl -big -ss -v --valley-seeking -p 0.05 --valley-depth 0.9 --dbkey mm10 --multi-test $f.tag.uniq.rgb.bed $f.tag.uniq.peak.sig.bed --out-boundary $f.tag.uniq.peak.sig.boundary.bed --out-half-PH $f.tag.uniq.peak.sig.halfPH.bed">>$sbatchname

echo "awk '{print \$1\"\t\"int((\$2+\$3)/2)-500\"\t\"int((\$2+\$3)/2)+500\"\t\"\$4\"\t\"\$5\"\t\"\$6}' $f.tag.uniq.peak.sig.halfPH.bed > $f.tag.uniq.peak.sig.halfPH.center.ext1k.bed">>$sbatchname
echo "perl ${ctkpath}/getMutationType.pl -t del $f.tag.uniq.mutation.txt $f.tag.uniq.del.bed">>$sbatchname
echo "perl ${ctkpath}/getMutationType.pl -t ins $f.tag.uniq.mutation.txt $f.tag.uniq.ins.bed">>$sbatchname
echo "perl ${ctkpath}/getMutationType.pl -t sub $f.tag.uniq.mutation.txt $f.tag.uniq.sub.bed">>$sbatchname
echo "perl ${ctkpath}/CIMS.pl -big -n 10 -p -outp $f.tag.uniq.del.pos.stat.txt -v -c $f.cache_del $f.tag.uniq.rgb.bed $f.tag.uniq.del.bed $f.tag.uniq.del.CIMS.txt">>$sbatchname
echo "perl ${ctkpath}/CIMS.pl -big -n 10 -p -outp $f.tag.uniq.ins.pos.stat.txt -v -c $f.cache_ins $f.tag.uniq.rgb.bed $f.tag.uniq.ins.bed $f.tag.uniq.ins.CIMS.txt">>$sbatchname
echo "perl ${ctkpath}/CIMS.pl -big -n 10 -p -outp $f.tag.uniq.sub.pos.stat.txt -v -c $f.cache_sub $f.tag.uniq.rgb.bed $f.tag.uniq.sub.bed $f.tag.uniq.sub.CIMS.txt">>$sbatchname

echo "awk '{if(\$9<=0.05) {print \$0}}' $f.tag.uniq.del.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > $f.tag.uniq.del.CIMS.s30.txt
cut -f 1-6 $f.tag.uniq.del.CIMS.s30.txt > $f.tag.uniq.del.CIMS.s30.bed">>$sbatchname
echo "awk '{if(\$9<=0.05) {print \$0}}' $f.tag.uniq.ins.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > $f.tag.uniq.ins.CIMS.s30.txt
cut -f 1-6 $f.tag.uniq.ins.CIMS.s30.txt > $f.tag.uniq.ins.CIMS.s30.bed">>$sbatchname
echo "awk '{if(\$9<=0.05) {print \$0}}' $f.tag.uniq.sub.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > $f.tag.uniq.sub.CIMS.s30.txt
cut -f 1-6 $f.tag.uniq.sub.CIMS.s30.txt > $f.tag.uniq.sub.CIMS.s30.bed">>$sbatchname
echo "awk '{print \$1\"\t\"\$2-10\"\t\"\$3+10\"\t\"\$4\"\t\"\$5\"\t\"\$6}' $f.tag.uniq.del.CIMS.s30.bed > $f.tag.uniq.del.CIMS.s30.21nt.bed" >>$sbatchname
echo "awk '{print \$1\"\t\"\$2-10\"\t\"\$3+10\"\t\"\$4\"\t\"\$5\"\t\"\$6}' $f.tag.uniq.ins.CIMS.s30.bed > $f.tag.uniq.ins.CIMS.s30.21nt.bed" >>$sbatchname
echo "awk '{print \$1\"\t\"\$2-10\"\t\"\$3+10\"\t\"\$4\"\t\"\$5\"\t\"\$6}' $f.tag.uniq.sub.CIMS.s30.bed > $f.tag.uniq.sub.CIMS.s30.21nt.bed" >>$sbatchname

echo "sbatch $sbatchname"
done
