#!/bin/bash
for i in {0..9}
do
mkdir Bjuipxml$i
echo "#!/bin/bash" > forscan$i.sh
echo "for file in \$(ls /media/tu/db/blastdb/Bju/v1.5/testfolds/pep$i)" >>forscan$i.sh
echo "do" >>forscan$i.sh
echo "if" >>forscan$i.sh 
echo "grep -q \$file /media/tu/db/blastdb/Bju/v1.5/testfolds/Bjuipxml$i/*.xml" >>forscan$i.sh
echo "then" >>forscan$i.sh
echo "continue" >>forscan$i.sh
#echo "$file exist"
echo "fi" >>forscan$i.sh
echo "./interproscan.sh -i /media/tu/db/blastdb/Bju/v1.5/testfolds/pep$i/\$file -f XML -d /media/tu/db/blastdb/Bju/v1.5/testfolds/Bjuipxml$i/ -dp -cpu 38 -iprlookup -goterms" >>forscan$i.sh
echo "done" >>forscan$i.sh
cp forscan$i.sh ./Bjuipxml$i
done
