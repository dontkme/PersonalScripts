#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2021
# Get Gencode translation sequences v1.01 2021/11/11
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
our $input;
# our $upstreml=5000;
# our $downstreml=5000;
# our $annot;
# our $sortID="ID";
# our $feature="mRNA";
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"i=s"=>\$input)
or die("[-]Error in command line arguments
  Usage: perl Getseqs [options] <-i string|input list file> <input FASTA file>
    options:
    [-o string|outprefix Default: getseqsOut]
    Note: Get Gencode translation sequences v1.01 2021/11/11.\n");

#############

open LIST, "< $input" or die ("[-] Error: Can't open input list file: $input.\n");
open OUT, "> $opfn.fa" or die ("[-] Error: Can't open or create $opfn.fa\n");

####### Start loading the sequences.
if (not $ARGV[0]) {
	die ("[-] Error: Not find a input Genome FASTA file.\n");
}
our $loadingstarttime=time();
# print @ARGV;
# if (@ARGV eq "") {
#     die ("[-] Error: Not find a input genome FASTA file.\n");
#     }
print "Start loading genomeic sequences.\n";

# our $seqi=0;
# our @Chrname=();
# our @Chrseq=();
our %Name2seq;
our $seqname="";
our $seqnum=0;

while(defined(our $seq = <>)){
#while(our $seq = <SEQFILENAME>){
	if ($seq =~ m/^.*>/) {
	$seq=~ m/^\s*>\s*(\S+)\s*(.*)/; #1.52 new RE based on BioPerl
	# say $1;
	my @seqnames=split('\|',$1);
	 $seqname= $seqnames[1]; # 1 is Gencode pc_translations.fa headline pattern. 
	#  say $seqname;
	$seqnum++;
	}else{
		$seq =~ s/\s//;
		# $seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
	 $Name2seq{$seqname}.=$seq;
	}
}

our $loadingendtime=time();
print "$seqnum Sequences\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;

#### Loading input list of Transcript IDs. fomat 1st name\t 2ed name\t gene name\n
print "Start loading input list.\n";
our $inputlineNum=0;
our $outputlineNum=0;
while(defined(our $listline=<LIST>)){
	$inputlineNum++;
	$listline =~ s/\R//g;
	my @tmplistline=split("\t", $listline);
	my $firstName=$tmplistline[0];
	my $secondName=$tmplistline[1];
	my $genename=$tmplistline[2];

	if(defined($firstName) && defined($secondName)){
		$outputlineNum++;
		print OUT ">$firstName:$genename\n";
		print OUT "$Name2seq{$firstName}\n";
		print OUT ">$secondName:$genename\n";
		print OUT "$Name2seq{$secondName}\n";
		
	}
}
print "Input list lines: $inputlineNum\nOutput: $outputlineNum\n";
close LIST;
close OUT;
