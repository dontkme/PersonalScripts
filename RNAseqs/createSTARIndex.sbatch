#!/bin/bash
#SBATCH --partition=bigmem2 		# partition nam
#SBATCH -c 10 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem=60G
#SBATCH --job-name STAR_index 		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /project2/xczhang/genomes/GENCODE_M25_STAR2.7_50_index --genomeFastaFiles /project2/xczhang/genomes/GENCODE_M24_GRCm38.p6/GRCm38.p6.genome.fa --sjdbGTFfile /project2/xczhang/genomes/GENCODE_M24_GRCm38.p6/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 50

