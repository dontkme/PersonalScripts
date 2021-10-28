#!/usr/bin/env perl
#AUTHORS
# Kaining Hu (c) 2021
# Get motif v0.100 2021/10/28
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
# use Parallel::ForkManager;
# our $MAX_processes=1;
# my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
# use Array::Utils qw(:all);
# use List::MoreUtils ':all';
# use Tie::Hash::Regex;

our $opfn="getMotifsOut";
my $verbose;
# our $upstreml=5000;
# our $downstreml=5000;
# our $annot = "";
# our $sortID = "gene_id";
# our $feature = "";
# our $inputlist;
our $Motif="";
# our $dj = 50;

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn, "f=s"=>\$Motif,"verbose"  => \$verbose)
or die("[-]Error in command line arguments
  Usage: perl GetMotif.pl [options] -f MotifPattern <input FASTA file(s)>
    options:
    [-o string|outprefix Default: getMotifsOut]
    [-f string|Motif RE. default: '']
    Note: Get Motif sequence name and Motif numbers v0.100 2021/10/28.\n");

if (not @ARGV) {
	die ("[-] Error: Not find a input FASTA file.\n");
}

if ($Motif eq ""){
    die ("[-] Error: Please enter the Motif pattern.\n");
}
open OUTPUTFILE , "> $opfn.$Motif.Out.txt" or die ("[-] Error: Can't open or create $opfn.$Motif.Out.txt\n");

our $loadingstarttime=time();
print "Start loading input sequence file(s): @ARGV \n";
our $Chri=0;
our $Chrname="";
our @Chrseq=();
our %Chrname2seq;

while(defined(our $seq = <>)){

	if ($seq =~ m/^.*>/) {

        $seq=~ m/^\s*>\s*(\S+)\s*(.*)/;
        # print "$1\n";
        $Chrname= $1;
        $Chri++;

	}else{
		$seq =~ s/\s//g;
		$seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;# snp replace
	    #$Chrseq[$Chri-1] .=$seq;
        $Chrname2seq{$Chrname} .=$seq;
		}
}
our $loadingendtime=time();
print "$Chri Sequence(s)\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;

our $starttime=time();

print "Motif: $Motif\n";
# open OUTPUTFILE , "> $opfn.$Motif.Out.txt" or die ("[-] Error: Can't open or create $opfn.$Motif.Out.txt\n");
print "Output file: $opfn.$Motif.Out.txt";
print "Running. Please wait for a minute.\n";
our $MotifSeqCount=0;
foreach my $key(%Chrname2seq){
    my $motifcount=0;
    my $tmpseq = $Chrname2seq{$key};
    while (defined $tmpseq && $tmpseq =~ m/$Motif/gc){

        $motifcount++;

    }
    if ($motifcount>0){
        $MotifSeqCount++;
        print OUTPUTFILE "$key\t$motifcount\n";
    }


}
our $endtime=time();
#say $starttime;
#say $endtime;
print "Total $MotifSeqCount sequence(s) have Motif($Motif).\n";
printf "All done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;



