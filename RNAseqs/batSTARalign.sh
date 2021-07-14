#!/bin/bash 
#DropseqD="/project2/xczhang/KnHu/biosoft/Drop-seq_tools-2.4.0" 
#STARindex="/project2/xczhang/genomes/mm10_STAR_index" 
STARindex="/project2/xczhang/genomes/mm10_STAR2.7.9_index"
 #metaFASTA="/project2/xczhang/genomes/mm10_meta/mm10/mm10.fasta"
#STARD="/software/STAR-2.6.1b-el7-x86_64/bin/STAR"
GTF="/project2/xczhang/genomes/mm10_meta/mm10/mm10.gtf"
#for file in $(ls *.unmapped.bam)
for file in $(ls *R1.fq.gz)
do
prefix=$(echo $file|sed 's/\(.*\)\R1\.fq\.gz/\1/')
#echo $prefix
f1="${prefix}R1.fq.gz"
f2="${prefix}R2.fq.gz"
#outputfolder="${prefix}_DropSeqRes"
#jobname="${prefix}.drop"
ncpu=4
jobname="${prefix}.STAR"
sbatchname="${prefix}.STAR2mm10.sbatch"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
echo "#SBATCH --partition=bigmem2" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=60G" >>$sbatchname

#echo "mkdir $outputfolder" >>$sbatchname
#echo "module load java" >>$sbatchname
#echo "sh ${DropseqD}/Drop-seq_alignment.sh -g $STARindex -r $metaFASTA -d $DropseqD -s $STARD -o $outputfolder $file " >>$sbatchname
echo "STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN $ncpu --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 3 --alignIntronMax 299999 --genomeDir $STARindex --sjdbGTFfile $GTF --outFileNamePrefix $prefix  --readFilesCommand zcat --readFilesIn $f1 $f2 " >>$sbatchname
echo "sbatch $sbatchname"
done
