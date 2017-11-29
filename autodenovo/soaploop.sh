#!/bin/bash
for((i=83;i<129;i=i+2)) 
do
mkdir $i
cp ex.config ./$i
cd $i
SOAPdenovo-127mer all -s ex.config -K $i -R -p 10  -k 45  -F -o B_K$i  1>ass.log 2>ass.err
cd ../
echo "$i done"
done
