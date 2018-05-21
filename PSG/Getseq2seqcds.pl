#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Combine 2 cds fasta by input files v1.0000 2018/05/17
# hukaining@gmail.com



use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
my $verbose;

our $cds1 = "/mnt/e/blast/db/BO/v1.0/B.oleracea_v1.0.scaffolds.cds";
our $cds2 = "/mnt/e/blast/db/Bnafra/BnafrCDS.fa";

our $cdsname="";
our %seqfa;
our $opdir = "./Getseq2seqcdsOut";


GetOptions("o=s" => \$opdir,"verbose"=>\$verbose,"c1=s"=>\$cds1,"c2=s"=>\$cds2)
or die("Error in command line arguments\nUsage: 	perl Getseq2seqcds.pl [options] -c1 <First_cds_fasta_file> -c2 <Second_cds_fasta_file> <Input 2seq name file>\n
	 options:\n
	 [-o string|output dir default: ./Getseq2seqcdsOut]\n
	 [-c1 string|CDS fasta file1]\n
	 [-c2 string|CDS fasta file2]\n
	 Note: Combine 2 cds fasta by input files v1.0000 2018/05/17\n");


our $loadingstarttime=time();

mkdir $opdir unless -d $opdir;

if(!open CDS1,"< $cds1"){
 die "$cds1 File not Found\n";
}


while(my $seq1 = <CDS1>){
  if ($seq1 =~ m/^.*>/) {
	$seq1=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
	 $cdsname = $1;
	$seqfa{$cdsname}=">$cdsname\n";
  }else {
    $seqfa{$cdsname}.=$seq1;
  }
}

close CDS1;

if(!open CDS2,"< $cds2"){
 die "$cds1 File not Found\n";
}

 $cdsname="";

while(my $seq1 = <CDS2>){
  if ($seq1 =~ m/^.*>/) {
	$seq1=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
	 $cdsname = $1;
	$seqfa{$cdsname}=$seq1;
  }else {
    $seqfa{$cdsname}.=$seq1;
  }
}

close CDS2;

print "Finished loading!\n";



our $count1=0;

while(defined(our $seq = <>)){
#while(our $seq = <SEQFILENAME>){
	my @sep = split/\t/,$seq;
	$sep[1] =~ s/\s//g;
	my $seqname="$sep[0]_2_"."$sep[1]";
	open OUTFA, "> $opdir/$seqname.fa" or die ("[-] Error: Can't open or creat $opdir/$seqname.fa\n");
	
  print OUTFA "$seqfa{$sep[0]}";
  print OUTFA "$seqfa{$sep[1]}";
  
	close OUTFA;
  $count1 ++;
			}
print "Finished combining 2CDS! $count1 fas in $opdir\n";
our $loadingendtime=time();
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;

