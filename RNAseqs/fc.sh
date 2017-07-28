#!/bin/bash
ls *.sort.bam >allfile.txt
echo "featureCounts -T 4 -t CDS -g gene_id -a Tair10.33.ZZQT.gtf -o ZZQ_fc_gene_CDS.txt" $(<allfile.txt)



