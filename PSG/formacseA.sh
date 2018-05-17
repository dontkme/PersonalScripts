#!/bin/bash
for file in $(ls /media/tu/db/HKN/PSG/A2BnaApaired/*.fa)
do
sfile=$(echo $file| sed 's/.*\/\(.*\)\.fa/\1/')
echo $sfile
if
grep -q $sfile /media/tu/db/HKN/PSG/A2BnaApaired/*NT.fasta
then
continue
fi
java -jar ./macse_v1.2.jar -prog alignSequences -seq $file 
done
