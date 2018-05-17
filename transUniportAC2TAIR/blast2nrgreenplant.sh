#!/bin/bash
time blastp -db /media/tu/db/nr/greenplant -task blastp-fast -outfmt 5 -query Senvy_insert.chr.modified_id.pep.fa -evalue 1e-3 -max_target_seqs 10 -out Bjuv1.5blastp2nr_greenplant_1e-3_max_10.xml -num_threads 39
