#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Convert 2 cds macse to KAKS calculator axt format v1.0000 2018/05/21
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
#use Time::HiRes 'time';
my $verbose;
our $axtname="";
our $seqfa="";
our $outaln = "./Out.axt";
if ($ARGV[0]) {
  $outaln="$ARGV[0].axt";
}

GetOptions("o=s" => \$outaln,"verbose"=>\$verbose)
or die("Error in command line arguments\nUsage: 	perl Macse2NT2axt.pl [options] <Input 2seq Macse out file>\n
	 options:\n
	 [-o string|output axt default: ./<Input 2seq Macse out file>.axt]\n
	 
	 Note: Convert 2 cds macse to KAKS calculator axt format v1.0000 2018/05/21\n");

#our $loadingstarttime=time();


while(defined(our $seq = <>)){
  if ($seq =~ m/^.*>/) {
	$seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
	 $axtname .= $1;
	
  }else {
    $seqfa.=$seq;
  }

}
  open OUTFA, "> $outaln" or die ("[-] Error: Can't open or creat $outaln\n");

  print OUTFA "$axtname\n";
  print OUTFA "$seqfa\n";
  
  close OUTFA;
  $axtname="";
  $seqfa="";