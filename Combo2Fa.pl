#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2021
# Transformat NMD combined result to Fasta file. (Combo2Fa)v0.15 2021/06/03
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
our $opfn="NMDSEFas";
my $verbose;

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
or die("[-]Error in command line arguments
  Usage: perl Combo2Fa [options] <input combined file>
    options:
	 [-o string|output prefix Default: Combined_fasta]
	
    Note: Transformat NMD combined result to Fasta file. (Combo2Fa)v0.15 2021/06/03\n");

############### Main ###########
open COMBOFA,"> $opfn.fa" or die ("[-] Error: Can't open or create $opfn.fa\n");
open COMBOFASIMPLE,"> $opfn.simple.fa" or die ("[-] Error: Can't open or create $opfn.simple.fa\n");

our $count = 1;
our @exonid = ();

while(defined(our $seq = <>)){

    if ($seq =~ m/^#/){
        next;

    }else{
        my @tmp = split("\t",$seq);
        my $originName = "$tmp[0].$tmp[1]:$tmp[2].$tmp[3].$count:$tmp[9]:$tmp[10]:$tmp[28]:$tmp[31]:$tmp[34]";
        chomp $originName;
        my $simpleexonid = "$tmp[0]:$tmp[3]:$count:$tmp[9]:$tmp[10]:$tmp[28]:$tmp[31]:$tmp[34]"; # [3] transcript_id [34] SEupstreamAApos [9] SE/US position type. [10] SE length. [28] ORF type. [31] NMD flag
        chomp $simpleexonid;
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

            my $SEtype="rm_or_add_SE";

        unless($tmp[21] eq "" or $tmp[22] eq ""){
            unless($orgininAASeqSimple =~/_/){
                
                if ($tmp[32] eq "USSEDS"){
                    $SEtype = "rm_SE";

                }elsif($tmp[32] eq "USDS"){
                    $SEtype = "Add_SE"
                }
                print COMBOFA ">$originName.original\n$originAASeq\n>$originName.$SEtype\n$afterAASeq\n";
                # print COMBOFASIMPLE ">$originNameSimple.original\n$orgininAASeqSimple\n>$originNameSimple.rm_or_add_SE\n$afterAASeqSimple\n";
                print COMBOFASIMPLE ">$simpleexonid.original\n$orgininAASeqSimple\n>$simpleexonid.$SEtype\n$afterAASeqSimple\n";
                $count++;
            }
        }


    }

	
}
close COMBOFA;
close COMBOFASIMPLE;
$count = $count -1;

say "Done. $count records. $opfn.fa $opfn.simple.fa";

