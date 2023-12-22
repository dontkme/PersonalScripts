#!/bin/bash
ncpu=4
mem=40G
#ref=/project2/xczhang/genomes/GENCODE_M25_IRFinder_index
#exe=/project2/xczhang/KnHu/biosoft/IRFinder-1.3.0/bin/IRFinder
ref=/project2/xczhang/genomes/GENCODE38_IRfinder1.31Index
#exe=/project2/xczhang/KnHu/biosoft/IRFinder-1.3.0/bin/IRFinder
exe=/project2/xczhang/KnHu/biosoft/IRFinder-1.3.1/bin/IRFinder
for i in $(ls -d *rep1.*.unsort.bam)
do
#echo $i
#echo "$ncpu $mem"
Prefix=$(echo $i|sed  's/\(.*\)\.rep1.*/\1/g' )
Suffix=$(echo $i|sed 's/.*rep1\.\(.*\.unsort.bam\)/\1/')
#echo $Prefix
file1=${Prefix}.rep1.$Suffix
file2=${Prefix}.rep2.$Suffix
#file2=XZ-XR-97S-${Prefix}_R2_001.fastq.gz
#echo $file1
#short=$(echo $Prefix|sed 's/\(.*\)_S.*/\1/g')
short=$Prefix.Pool
#echo $short

echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short.IR2GCh38" >>$short.sbatch
echo "#SBATCH --output=$short.out" >>$short.sbatch
echo "#SBATCH --error=$short.err" >>$short.sbatch
#echo "#SBATCH --partition=bigmem2" >>$short.sbatch
echo "#SBATCH --account=pi-xczhang" >>$short.sbatch
echo "#SBATCH --partition=caslake" >>$short.sbatch
echo "#SBATCH --ntasks=${ncpu}" >>$short.sbatch
echo "#SBATCH --mem=${mem}" >>$short.sbatch

echo "module load gcc">>$short.sbatch
#echo "$exe -t $ncpu -r $ref -a 'none' -d $short.od $file1 $file2">>$short.sbatch
echo "$exe -t $ncpu -m BAM -r $ref -a 'none' -d $short.od <(samtools cat -@ $ncpu $file1 $file2)">>$short.sbatch
echo "sbatch $short.sbatch"

done
