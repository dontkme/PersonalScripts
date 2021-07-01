#!/bin/bash
DropseqD="/project2/xczhang/KnHu/biosoft/Drop-seq_tools-2.4.0"
#STARindex="/project2/xczhang/genomes/mm10_STAR_index"
STARindex="/project2/xczhang/genomes/mm10_STAR2.7_index"
metaFASTA="/project2/xczhang/genomes/mm10_meta/mm10/mm10.fasta"
STARD="/software/STAR-2.6.1b-el7-x86_64/bin/STAR"
for file in $(ls *.unmapped.bam)
do
prefix=$(echo $file|sed 's/\(.*\)\.unmapped\.bam/\1/')
#echo $prefix
outputfolder="${prefix}_DropSeqRes"
jobname="${prefix}.drop"
sbatchname="${prefix}.DropsAlig.sbatch"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
echo "#SBATCH --partition=bigmem2" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=1" >>$sbatchname
echo "#SBATCH --mem=60G" >>$sbatchname

echo "mkdir $outputfolder" >>$sbatchname
echo "module load java" >>$sbatchname
#echo "sh ${DropseqD}/Drop-seq_alignment.sh -g $STARindex -r $metaFASTA -d $DropseqD -s $STARD -o $outputfolder $file " >>$sbatchname
echo "sh ${DropseqD}/Drop-seq_alignment.sh -g $STARindex -r $metaFASTA -d $DropseqD -o $outputfolder $file " >>$sbatchname
echo "sbatch $sbatchname"
done
