#!/bin/bash
#SBATCH --job-name=fomatnr
#SBATCH --output=formatnr.out
#SBATCH --error=formatnr.err

#SBATCH --partition=broadwl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G

time ~/scratch-midway2/I-TASSER5.1/blast/bin/formatdb -t nr_old -i nr -p T -v 10000
