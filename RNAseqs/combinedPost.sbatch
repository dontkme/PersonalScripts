#!/bin/bash
#SBATCH --job-name=combinedPost
#SBATCH --output=combinedPost.out
#SBATCH --error=combinedPost.err
#SBATCH --partition=bigmem2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
source activate /home/kaininghu/biosoft/rmats_turbo_v4_1_1/conda_envs/rmats
python3 /home/kaininghu/xczhang/KnHu/biosoft/rmats-turbo-kutscherae-individual-counts/rmats.py --b1 AllsampleBAMs.txt --od ./Allsamples.od --tmp ./Allsamples.tmp  --anchorLength 1 --readLength 101 --variable-read-length --gtf /project2/xczhang/genomes/GENCODE_M24_GRCm38.
p6/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf -t paired --task post --nthread 4 --statoff
#python3 /home/kaininghu/xczhang/KnHu/biosoft/rmats-turbo-kutscherae-individual-counts/rmats.py --b1 AllsampleBAMs.txt --od ./Allsamples.od --tmp ./Allsamples.tmp  --anchorLength 1 --readLength 101 --variable-read-length --gtf /project2/xczhang/genomes/GENCODE_M24_GRCm38.p6/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf -t paired --task post --nthread 4 --statoff
