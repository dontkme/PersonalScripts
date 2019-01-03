#!/bin/bash
#bgzip -c ZHYA02_5380005_5387257_addheader.vcf >ZHYA02_5380005_5387257_addheader.vcf.gz
#tabix -p vcf ZHYA02_5380005_5387257_addheader.vcf.gz
for line in $(cat sampleIDs505.txt)
do
    echo $line
    samtools faidx /media/db/Bnafr/Brassica_napus_v4.1.chromosomes.fa chrA02:5380005-5387257 |vcf-consensus ZHYA02_5380005_5387257_addheader.vcf.gz -s $line >out506_$line.fa
done

