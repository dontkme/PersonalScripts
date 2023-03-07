#!/bin/bash
ncpu=4
mem=60G
ref=/project2/xczhang/genomes/GENCODE_M25_IRFinder_index
exe=/project2/xczhang/KnHu/biosoft/IRFinder-1.3.0/bin/IRFinder
for i in $(ls *.unsort.bam)
do
#echo $i
#echo "$ncpu $mem"
Prefix=$(echo $i|sed  's/\(.*\)\.unsort\.bam/\1/g' )
#echo $Prefix
#file1=XZ-XR-97S-${Prefix}_R1_001.fastq.gz
#file2=XZ-XR-97S-${Prefix}_R2_001.fastq.gz
file1=${Prefix}.unsort.bam
#echo $file1
#short=$(echo $Prefix|sed 's/\(.*\)_S.*/\1/g')
short=$Prefix
#echo $short

echo "#!/bin/bash" >$short.sbatch
echo "#SBATCH --job-name=$short.IR2GCM25" >>$short.sbatch
echo "#SBATCH --output=$short.out" >>$short.sbatch
echo "#SBATCH --error=$short.err" >>$short.sbatch
echo "#SBATCH --partition=bigmem2" >>$short.sbatch
echo "#SBATCH --ntasks=${ncpu}" >>$short.sbatch
echo "#SBATCH --mem=${mem}" >>$short.sbatch

echo "module load gcc">>$short.sbatch
#echo "$exe -t $ncpu -r $ref -a 'none' -d $short.od $file1 $file2">>$short.sbatch
echo "$exe -t $ncpu -m BAM -r $ref -a 'none' -d $short.od $file1">>$short.sbatch
echo "sbatch $short.sbatch"

done
