awk '{if($6==0&&$7==0){print $0}}' Bnafr.denovo.cds.status|awk '{for(i=2;i<=NF;i++){if($i>0){a[$1]=i-1}}}END{for(i in a){print i,a[i]}}' >Bnafr.new.orphan.gene.candidate
