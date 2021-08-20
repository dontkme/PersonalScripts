#!/bin/bash
ncpu=8
for file in $(ls *.bam)
do 
echo $file >$file.sample.txt
#outdir="./sb_${file}.od"
outdir="./sb_combined.od"
gtf="/project2/xczhang/genomes/GENCODE_M24_GRCm38.p6/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf"
prefix=$file
jobname="${prefix}.rmats"
sbatchname="${prefix}.pairedrmatsprep.sbatch"
echo "#!/bin/bash" >$sbatchname
echo "#SBATCH --job-name=$jobname" >>$sbatchname
echo "#SBATCH --output=$jobname.out" >>$sbatchname
echo "#SBATCH --error=$jobname.err" >>$sbatchname
echo "#SBATCH --partition=bigmem2" >>$sbatchname
echo "#SBATCH --ntasks=1" >>$sbatchname
echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
echo "#SBATCH --mem=60G" >>$sbatchname

echo "source activate /home/kaininghu/biosoft/rmats_turbo_v4_1_1/conda_envs/rmats" >>$sbatchname
#echo "python3 /home/kaininghu/xczhang/KnHu/biosoft/rmats-turbo-kutscherae-individual-counts/rmats.py --b1 $file.sample.txt --od ${outdir} --tmp ${outdir}/.tmp --anchorLength 1 --readLength 100 --variable-read-length --gtf ${gtf} -t single --task prep --nthread 8 --statoff" >>$sbatchname
echo "python3 /home/kaininghu/xczhang/KnHu/biosoft/rmats-turbo-kutscherae-individual-counts/rmats.py --b1 $file.sample.txt --od ${outdir} --tmp ${outdir}/.tmp --anchorLength 1 --readLength 100 --variable-read-length --gtf ${gtf} -t paired --task prep --nthread 8" >>$sbatchname
echo "sbatch $sbatchname"
done
