#!/bin/bash
#for i in Bnafr Ath11 BniB Bol_HDEM Brapa3.0 Tpa
for i in BolHDEM 
do blat ../${i}.dna.fa ../Bnafr.cds.fa Bnafr.${i}.dna.psl
done
