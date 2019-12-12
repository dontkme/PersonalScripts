#!/bin/bash
vcfFile=$(ls *allindel.vcf.gz)
echo "#!/bin/bash" >indelcmd.sh
for ((i=1;i<11;++i))
do
#echo $i
ii=$(printf "%02d" $i) 
#echo $ii
#echo "grep -P \"^A${ii}\t\"  $vcfFile> chrA$ii.indel.vcf" >> indelcmd.sh
echo "bcftools view $vcfFile -r scaffoldA$ii -o A$ii.indel.vcf" >> indelcmd.sh
#echo "grep -P \"^chrA${ii}_random\"  $vcfFile> chrA${ii}_random.indel.vcf" >> indelcmd.sh
done
for ((i=1;i<10;++i))
do
#echo $i
ii=$(printf "%02d" $i)
#echo $ii
#echo "grep -P \"^C$ii\t\"  $vcfFile> chrC$ii.indel.vcf" >> indelcmd.sh
#echo "grep -P \"^chrC${ii}_random\"  $vcfFile> chrC${ii}_random.indel.vcf" >> indelcmd.sh
echo "bcftools view $vcfFile -r scaffoldC$ii -o C$ii.indel.vcf" >> indelcmd.sh
done
#echo "grep -P \"^chrAnn_random\"  $vcfFile> chrAnn_random.indel.vcf" >> indelcmd.sh
#echo "grep -P \"^chrCnn_random\"  $vcfFile> chrCnn_random.indel.vcf" >> indelcmd.sh
#echo "grep -P \"^chrUnn_random\"  $vcfFile> chrUnn_random.indel.vcf" >> indelcmd.sh
#echo "grep -P \"^scaffold\"  $vcfFile> scaffold.indel.vcf" >> indelcmd.sh
