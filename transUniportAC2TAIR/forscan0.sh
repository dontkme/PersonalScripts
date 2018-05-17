#!/bin/bash
for file in $(ls /media/tu/db/blastdb/Bju/v1.5/testfolds/pep0)
do
if
grep -q $file /media/tu/db/blastdb/Bju/v1.5/testfolds/Bjuipxml0/*.xml
then
continue
fi
./interproscan.sh -i /media/tu/db/blastdb/Bju/v1.5/testfolds/pep0/$file -f XML -d /media/tu/db/blastdb/Bju/v1.5/testfolds/Bjuipxml0/ -dp -cpu 38 -iprlookup -goterms
done
