#!/bin/bash
vcfFile=AS2127all_AS94570allBY.allsnp.vcf
echo "#!/bin/bash" >snpcmd.sh
for ((i=1;i<11;++i))
do
#echo $i
ii=$(printf "%02d" $i) 
#echo $ii
echo "grep -P \"^A${ii}\t\"  $vcfFile> chrA$ii.snp.vcf" >> snpcmd.sh
#echo "grep -P \"^chrA${ii}_random\"  $vcfFile> chrA${ii}_random.snp.vcf" >> snpcmd.sh
done
for ((i=1;i<10;++i))
do
#echo $i
ii=$(printf "%02d" $i)
#echo $ii
echo "grep -P \"^C$ii\t\"  $vcfFile> chrC$ii.snp.vcf" >> snpcmd.sh
#echo "grep -P \"^chrC${ii}_random\"  $vcfFile> chrC${ii}_random.snp.vcf" >> snpcmd.sh
done
#echo "grep -P \"^chrAnn_random\"  $vcfFile> chrAnn_random.snp.vcf" >> snpcmd.sh
#echo "grep -P \"^chrCnn_random\"  $vcfFile> chrCnn_random.snp.vcf" >> snpcmd.sh
#echo "grep -P \"^chrUnn_random\"  $vcfFile> chrUnn_random.snp.vcf" >> snpcmd.sh
echo "grep -P \"^scaffold\"  $vcfFile> scaffold.snp.vcf" >> snpcmd.sh
