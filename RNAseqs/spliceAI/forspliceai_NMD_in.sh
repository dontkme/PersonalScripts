#!/bin/bash
spliceD=4000
ref=GRCh38.p13.genome.fa
spliceA=grch38
for file in $(ls Test3_*NMD_in*chr*.vcf)
do

#echo $file
chr=$(echo $file|sed  's/.*\(chr.*\)\.vcf/\1/g')
#echo $chr

outfile="Test3_humanunion_fixbug_rmNMD_GTF.NMD_inAS449.intron_region2.gnomad.${chr}.spliceai.D${spliceD}.vcf"
#echo $outfile
echo "#!/bin/bash" >${file}.sbatch
echo "#SBATCH --job-name=NMD_in_spliceai.${chr}.D${spliceD}">>${file}.sbatch
echo "#SBATCH --output=NMD_in_spliceai.${chr}.D${spliceD}.out">>${file}.sbatch
echo "#SBATCH --error=NMD_in_spliceai.${chr}.D${spliceD}.err">>${file}.sbatch
echo "#SBATCH --nodes=1">>${file}.sbatch
echo "#SBATCH --partition=gpu2">>${file}.sbatch
echo "#SBATCH --ntasks=1">>${file}.sbatch
echo "#SBATCH --gres=gpu:1">>${file}.sbatch
echo "#SBATCH -w midway2-0664">>${file}.sbatch

echo "source activate /project2/xczhang/KnHu/condaenv/spliceai">>${file}.sbatch
echo "module load cuda/11.2">>${file}.sbatch
echo "spliceai -I $file  -O $outfile  -A $spliceA -R $ref -D $spliceD">>${file}.sbatch


done
