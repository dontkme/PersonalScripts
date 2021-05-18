#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2018
# Get Sequences from GFF with genome v1.1001 2018/09/12
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="getseqsOut";
my $verbose;
our $upstreml=5000;
our $downstreml=5000;
our $annot;
our $sortID="ID";
our $feature="mRNA";
our $inputlist="";

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"u=i"=>\$upstreml,"d=i"=>\$downstreml, "g=s"=>\$annot,"f=s"=>\$feature,"s=s"=>\$sortID, "i=s"=>\$inputlist)
or die("[-]Error in command line arguments
  Usage: perl Get3UTR [options] <-g string|GTF annoation file> <input all 3 UTR FASTA file(s)>
    options:
    [-o string|outprefix Default: getseqsOut]
    [-s string|Specify attribute type in GFF annotation for sorting. default: gene_id]
    [-f string|Specify feature type in GFF annotation.default: CDS]
    [-u int|upstream length Default: 5000]
    [-d int|downstream length Default: 5000]
    [-i input list]
	 
    Note: Get UTR Sequences  v1.0000 2021/05/17.\n");

###################sub TRseq##########
 


sub TRseq($)
{
	my ($pinseq) = @_;
	#say $pinseq;
	my $pinseqtr = reverse $pinseq;
	#say $pinseqtr;
	 $pinseqtr =~ tr/ACGTacgt/TGCAtgca/;
	 #say $pinseqtr;
	 return  $pinseqtr;
}
##################TRseq End#############
if (not $ARGV[0]) {
	die ("[-] Error: Not find a input Genome FASTA file.\n");
}
our $loadingstarttime=time();
# print @ARGV;
# if (@ARGV eq "") {
#     die ("[-] Error: Not find a input genome FASTA file.\n");
#     }
print "Start loading genomeic sequences.\n";

our $Chri=0;
our @Chrname=();
our @Chrseq=();
our %Chrid2seq;
#@ARGV = qw#''  Not_Find_a_File#;
#say @ARGV;
#say $0;
while(defined(our $seq = <>)){
    # chomp $seq;

	if ($seq =~ m/^.*>/) {
        $seq =~ m/^\s*>\s*((\S+)\s*(.*))\n/;
        print "$3\n";
        $Chrname[$Chri]= $3;
        $Chri++;
	}else{
		$seq =~ s/\s//;
		$seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
	    #$Chrseq[$Chri-1] .=$seq;
        $Chrid2seq{$Chrname[$Chri-1]} .=$seq;
	}
}

    # for (our $i=0;$i<$Chri-1;$i++){
    #     $Chrid2seq{$Chrname[$i]}=$Chrseq[$i];
    # }

#close SEQFILENAME;
our $loadingendtime=time();
print "$Chri Sequences\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;
#  sleep 2000;

open INPUTLIST, "< $inputlist" or die ("[-] Error: Can't open list file: $inputlist.\n");
our %key2seq;

while (our $inputline = <INPUTLIST>){
    chomp $inputline;
    # print "We have matched " if /$inputline/ ~~ %Chrid2seq; 
    foreach our $tmpkey (keys %Chrid2seq){
        if ($tmpkey =~ m/$inputline/){
            # say $tmpkey;
            if (defined $Chrid2seq{$tmpkey}){
                my $tmpseq = $Chrid2seq{$tmpkey};
                # say $tmpseq;
                $key2seq{$inputline} = $tmpseq;
            }
        }
    }
}
close INPUTLIST;

open INPUTLIST, "< $inputlist" or die ("[-] Error: Can't open list file: $inputlist.\n");

open OUT, "> $opfn.fa" or die ("[-] Error: Can't open or creat $opfn.fa\n");
while (my $inputline= <INPUTLIST>){
    chomp $inputline;
    # print "We have matched " if /$inputline/ ~~ %Chrid2seq; 
   
    # my $tmpkey = $key2seq{$inputline};
    # my $finalseq = $Chrid2seq{$tmpkey};
    # say $inputline;
    my $finalseq = $key2seq{$inputline};
    # say $finalseq;
    if(defined $finalseq){
    print OUT ">$inputline.3UTR\n";
    print OUT "$finalseq\n";
    }
}