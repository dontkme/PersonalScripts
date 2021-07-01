#!/bin/bash
#SBATCH --partition=bigmem2 		# partition nam
#SBATCH -c 2 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem=60G
#SBATCH --job-name STAR_index 		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir /project2/xczhang/genomes/mm10_STAR2.7_index --genomeFastaFiles /project2/xczhang/genomes/mm10_meta/mm10/mm10.fasta --sjdbGTFfile /project2/xczhang/genomes/mm10_meta/mm10/mm10.gtf --sjdbOverhang 99
