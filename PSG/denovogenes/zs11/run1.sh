#!/bin/bash
for i in Bnazs11 BniB BolHDEM Brapa3 Tpa Ath11
do blat ../${i}.cds.fa ../Bnazs11.cds.fa Bnazs11.${i}.cds.psl
done
