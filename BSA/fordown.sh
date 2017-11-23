#!/bin/bash
echo "#!/bin/bash">downlist.sh
while read i 
do
echo "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR320/$i/$i.sra">>downlist.sh
echo "echo "$i.sra done."">>downlist.sh
done < 'SRR_Acc_List.txt'
