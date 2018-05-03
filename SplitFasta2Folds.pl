#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# split fasta file to saparate files.to different folds (SpliteFasta2Folds)v1.0000 2018/05/03
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
our $opfn="";
my $verbose;
#our $upstreml=3000;
our $opdir='./';
our $suffix='';
our $perfoldseqnum=200000;
our $foldnum=0;
#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opdir,"verbose"=>\$verbose,"n=i" => \$perfoldseqnum,"s=s"=>\$suffix)
or die("Error in command line arguments\nUsage: perl splitFasta.pl [-o string |outfiledir][-s string |suffix][-n int |perfold seq number]<input_fasta>\n");


our $loadingstarttime=time();

print "Start loading fasta sequence(s).\n";
#if(!open SEQFILENAME,"< $seqfilename"){
# die "File not Found\n";
#}
#if(!open SEQFILENAME,"< $_"){
# die "File not Found\n";
#}
if ($perfoldseqnum<=0) {die ("-n perfold seq number should larger than zero.\n")};
our $Chri=0;
our @Chrname=();
#our @Chrseq=();
#@ARGV = qw#''  Not_Find_a_File#;
#say @ARGV;
#say $0;
while(defined(our $seq = <>)){
#while(our $seq = <SEQFILENAME>){
	if ($seq =~ m/^.*>/) {
	$seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
	close OUTFA;
	$foldnum =int($Chri/$perfoldseqnum);
	my $comdir="$opdir$foldnum"; 
	mkdir $comdir unless -d $comdir;
	my $seqname=$1;
	open OUTFA, "> $comdir/$seqname$suffix" or die ("[-] Error: Can't open or creat $comdir/$seqname$suffix\n");
  print OUTFA ">$1\n";
	
	 $Chrname[$Chri]= $1;
	$Chri++;
	}else{
		#$seq =~ s/\s//;
	 #$seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
	 #$Chrseq[$Chri-1] .=$seq;
	 print OUTFA "$seq";
			}
}


close OUTFA;
our $loadingendtime=time();
$foldnum++;
print "$Chri Sequences\n$perfoldseqnum perfold\n$foldnum Fold(s)\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;
