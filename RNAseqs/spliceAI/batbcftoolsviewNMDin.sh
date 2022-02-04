#!/bin/bash
for i in {1..22} X Y
do
echo $i 
bcftools view -R Test3_humanunion_fixbug_rmNMD_GTF.NMD_exAS968.intron_region2.chr.txt ../gnomad.genomes.v3.1.2.sites.chr${i}.vcf.bgz >Test3_humanunion_fixbug_rmNMD_GTF.NMD_exAS968.intron_region2.gnomad.chr${i}.vcf
done
