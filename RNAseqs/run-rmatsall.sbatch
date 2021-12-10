#!/bin/bash
#SBATCH --job-name=runrmatsall2
#SBATCH --output=runrmatsall2.out
#SBATCH --error=runrmatsall2.err
#SBATCH --partition=bigmem2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
source activate /home/kaininghu/biosoft/rmats_turbo_v4_1_1/conda_envs/rmats
export LD_LIBRARY_PATH=/home/kaininghu/biosoft/rmats_turbo_v4_1_1/conda_envs/rmats/lib:$LD_LIBRARY_PATH
python3 /home/kaininghu/xczhang/KnHu/biosoft/rmats-turbo-kutscherae-individual-counts/rmats.py --b1 Pos.txt --b2 Neg.txt --od ./AllsamplesR.od --tmp ./AllsamplesR.tmp --readLength 58 --variable-read-length --gtf /project2/xczhang/genomes/Mus_musculus_Ensembl_GRCm38/Ensembl/GRCm38/Annotation/Genes/genes.gtf -t paired --task both --nthread 4 
