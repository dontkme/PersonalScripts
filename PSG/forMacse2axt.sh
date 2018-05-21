#!/bin/bash
for file in $(ls *NT.fasta)
do
perl ../Macse2NT2axt.pl -o $file.1.axt $file
done
