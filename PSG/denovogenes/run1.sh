#!/bin/bash
for i in Bnafr Ath11 BniB BolHDEM Brapa3 Tpa
do blat ../${i}.cds.fa ../Bnafr.cds.fa Bnafr.${i}.cds.psl
done
