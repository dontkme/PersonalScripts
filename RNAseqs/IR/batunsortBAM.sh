#!/bin/bash
ncpu=4
for i in $(ls *.bam)
do	
prefix=$(echo $i|sed 's/\(.*\)\.bam/\1/g')
jobname="${prefix}.unsort"
sbatchname="${prefix}.unsort.sbatch"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
echo "#SBATCH --account=pi-xczhang" >>$sbatchname
echo "#SBATCH --partition=caslake" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=20G" >>$sbatchname
	#echo $i
	#prefix=$(echo $i|sed 's/\(.*\)\.bam/\1/g')
       #echo $prefix
echo "samtools sort -@ $ncpu -n -o $prefix.unsort.bam $i" >> $sbatchname
echo "sbatch $sbatchname"
done
