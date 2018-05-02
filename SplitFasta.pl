#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# split fasta file to saparate files. (EAHelitron)v1.0000 2018/05/02
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
our $suffix='.fasta';
#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opdir,"verbose"=>\$verbose,"s=s"=>\$suffix)
or die("Error in command line arguments\nUsage: perl splitFasta.pl [-o|outfiledir string][-s|suffix string] <input_fasta>\n");


our $loadingstarttime=time();

print "Start loading fasta sequence(s).\n";
#if(!open SEQFILENAME,"< $seqfilename"){
# die "File not Found\n";
#}
#if(!open SEQFILENAME,"< $_"){
# die "File not Found\n";
#}
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
	my $seqname=$1;
	open OUTFA, "> $opdir/$seqname$suffix" or die ("[-] Error: Can't open or creat $opdir/$seqname$suffix\n");
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
print "$Chri Sequences\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;