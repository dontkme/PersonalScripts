#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2021
# Transformat NMD combined result to Fasta file. (Combo2Fa)v0.10 2021/06/01
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
# use Time::HiRes 'time';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="EAHeli_out";
my $verbose;

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
or die("[-]Error in command line arguments
  Usage: perl Combo2Fa [options] <input combined file>
    options:
	 [-o string|output prefix Default: Combined_fasta]
	
    Note: Transformat NMD combined result to Fasta file. (Combo2Fa)v0.10 2021/06/01\n");

############### Main ###########
open COMBOFA,"> $opfn.fa" or die ("[-] Error: Can't open or create $opfn.fa\n");
open COMBOFASIMPLE,"> $opfn.simple.fa" or die ("[-] Error: Can't open or create $opfn.simple.fa\n");

our $count = 0;
our @exonid = ();

while(defined(our $seq = <>)){

    if ($seq =~ m/^#/){
        next;

    }else{
        my @tmp = split("\t",$seq);
        my $originName = "$tmp[0].$tmp[1].$tmp[2].$tmp[3].$count";
        my $simpleexonid = "$tmp[0]:$tmp[3].$count"; # [3] transcript_id
        my $originNameSimple = $originName;
            $originNameSimple =~ s/@//g;
        my $originAASeq = $tmp[21];
        my $orgininAASeqSimple = $originAASeq;
            $orgininAASeqSimple =~ s/_+$//;
            # $orgininAASeqSimple = $orgininAASeqSimple;
            $orgininAASeqSimple =~ s/\d+$//;

        my $afterAASeq = $tmp[22];
        my $afterAASeqSimple = $afterAASeq;
            $afterAASeqSimple =~ s/_.*//;
            $afterAASeqSimple =~ s/\d*$//;

        unless($tmp[21] eq "" or $tmp[22] eq ""){
            unless($orgininAASeqSimple =~/_/){
                print COMBOFA ">$originName.original\n$originAASeq\n>$originName.rm_or_add_SE\n$afterAASeq\n";
                # print COMBOFASIMPLE ">$originNameSimple.original\n$orgininAASeqSimple\n>$originNameSimple.rm_or_add_SE\n$afterAASeqSimple\n";
                print COMBOFASIMPLE ">$simpleexonid.original\n$orgininAASeqSimple\n>$simpleexonid.rm_or_add_SE\n$afterAASeqSimple\n";
                $count++;
            }
        }


    }

	
}
close COMBOFA;
close COMBOFASIMPLE;

say "Done. $count records. $opfn.fa $opfn.simple.fa";

