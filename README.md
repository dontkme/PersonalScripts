# PersonalScripts
My PersonalScripts.  

transBnafrat2atpro.pl: Used for blast results(outfmt6) add descriptions.

batentry: blastdbcmd batch_entry perl script.

forsepSTRUCTURE.sh  bat STRUCTURE cmd maker.

RandomDNA_rate.pl 

Creat a random DNA seq with your ATGC rate v1.0000 2018/03/08

    Usage:  perl RandomDNA.pl [options]

         options:

         [-o string|output prefix default: rnddna100]

         [-n string|Add fasta header default: "randomDNA" ]

         [-l int|Length of the DNA. default: 100]

         [-r string|A,T,G,C rate, sep by comma. default:1,1,1,1]
         
         
         
CountATGC.pl 

Count Seq ATGC (CountATGC)v1.1000 2018/03/09
    
    Usage: perl countATGC [-o outfileprefix] <inputfile>
