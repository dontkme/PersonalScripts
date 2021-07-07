#!/bin/bash
MG=200
DropseqPATH=/project2/xczhang/KnHu/biosoft/Drop-seq_tools-2.4.0
for i in $(ls -d *Res)
do
short=${i}.${MG} 
echo "#!/bin/bash" > $short.sbatch
echo "#SBATCH --job-name=$short.MG.$MG" >>$short.sbatch
echo "#SBATCH --output=$short.out" >>$short.sbatch
echo "#SBATCH --error=$short.err" >>$short.sbatch
echo "#SBATCH --partition=broadwl" >>$short.sbatch
echo "#SBATCH --ntasks=1" >>$short.sbatch
echo "#SBATCH --mem=8G" >>$short.sbatch
echo "module load java" >>$short.sbatch
echo "cd $i" >>$short.sbatch
#echo "java -jar -Xmx8G ~/xczhang/KnHu/biosoft/picard. -F1 $file1 -F2 $file2 -O $short.unmapped.bam -SM $short -RG $short" >>$short.sbatch
#echo "${DropseqPATH}/DigitalExpression -I final.bam -O DGE.MG${MG}.out.txt -SUMMARY DGE.MG${MG}.summary.txt -MIN_NUM_GENES_PER_CELL ${MG}" >>$short.sbatch
echo "${DropseqPATH}/DigitalExpression I=final.bam O=DGE.MG${MG}.out.txt SUMMARY=DGE.MG${MG}.summary.txt MIN_NUM_GENES_PER_CELL=${MG}" >>$short.sbatch
echo "sbatch $short.sbatch"
done
