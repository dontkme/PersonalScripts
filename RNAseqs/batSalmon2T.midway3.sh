#!/bin/bash 
#DropseqD="/project2/xczhang/KnHu/biosoft/Drop-seq_tools-2.4.0" 
#STARindex="/project2/xczhang/genomes/mm10_STAR_index" 
#STARindex="/project2/xczhang/genomes/GENCODE_M25_STAR2.7_index"
#STARindex="/project2/xczhang/genomes/GENCODE38_GRCh38.p13_ALL_STAR2.7.9_index"
SMindex="/project/xczhang/Genomes/GENCODE43_salmonTrans_Index/GENCODE43_Trans_salmon_index"
#metaFASTA="/project2/xczhang/genomes/mm10_meta/mm10/mm10.fasta"
SMPATH="/home/kaininghu/biosoft/salmon-latest_linux_x86_64/bin/"
#STARD="/software/STAR-2.6.1b-el7-x86_64/bin/STAR"
#GTF="/project2/xczhang/genomes/GENCODE_M24_GRCm38.p6/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf"
GTF="/project2/xczhang/genomes/GENCODE_38_GRCh38.p13/gencode.v38.primary_assembly.annotation.gtf"
#for file in $(ls *.unmapped.bam)
for file in $(ls *R1.clean.fq.gz)
do
prefix=$(echo $file|sed 's/\(.*\)\.R1\.clean\.fq\.gz/\1/')
#echo $prefix
f1="${prefix}.R1.clean.fq.gz"
f2="${prefix}.R2.clean.fq.gz"
prefix2=$prefix.GCh38
prefix2=$prefix.GCh38
#outputfolder="${prefix}_DropSeqRes"
#jobname="${prefix}.drop"
ncpu=10
jobname="${prefix2}.SM2T"
sbatchname="${prefix2}.SM2T.sbatch"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
echo "#SBATCH --account=pi-xczhang" >>$sbatchname
echo "#SBATCH --partition=caslake" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=20G" >>$sbatchname

#echo "mkdir $outputfolder" >>$sbatchname
#echo "module load java" >>$sbatchname
#echo "sh ${DropseqD}/Drop-seq_alignment.sh -g $STARindex -r $metaFASTA -d $DropseqD -s $STARD -o $outputfolder $file " >>$sbatchname
#echo "STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN $ncpu --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 3 --alignIntronMax 299999 --genomeDir $STARindex --sjdbGTFfile $GTF --outFileNamePrefix $prefix  --readFilesCommand zcat --readFilesIn $f1 $f2 " >>$sbatchname
#echo "STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN $ncpu --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 3 --alignIntronMax 299999 --genomeDir $STARindex --sjdbGTFfile $GTF --outFileNamePrefix $prefix2  --readFilesCommand zcat --readFilesIn $f1 $f2 " >>$sbatchname
#echo "samtools index -@ $ncpu ${prefix2}Aligned.sortedByCoord.out.bam" >>$sbatchname
echo "$SMPATH/salmon quant -i $SMindex -l A -1 ${f1} -2 ${f2} -p $ncpu --validateMappings -o Tquants/${prefix}_quant" >>$sbatchname
echo "sbatch $sbatchname"
done
