#!/bin/bash
for i in $(ls /Users/kaininghu/HKN/0Postdoc/202201/IR/IRFinder_res/*rep1.od/*-dir.txt)
do
samplename=$(echo $i|sed 's/.*\(S0.*\)\.rep1.*/\1/')
#echo $samplename
cmpname=WT_$samplename
echo $cmpname
mkdir $cmpname
cp ./experiment.txt ./$cmpname
cp ./filePaths.txt ./$cmpname
file1=$i
file2=$(echo $i|sed 's/rep1/rep2/')
samplenamerep1=$samplename.rep1
samplenamerep2=$samplename.rep2
echo "$file1
$file2" >>./$cmpname/filePaths.txt
echo "$samplenamerep1	KO
$samplenamerep2	KO" >>./$cmpname/experiment.txt

cd $cmpname
Rscript ../AutoIR.R
cd ..

done
