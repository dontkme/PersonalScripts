#!/bin/bash
for i in Bnazs11 Ath11 BniB BolHDEM Brapa3 Tpa
#for i in BolHDEM 
do blat ../${i}.dna.fa ../Bnazs11.cds.fa Bnazs11.${i}.dna.psl
done
