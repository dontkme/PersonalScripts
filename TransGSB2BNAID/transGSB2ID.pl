#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Trans GSBID to BNAID. v1.0001 2018/07/21
# hukaining@gmail.com



use strict;
use warnings;
use 5.0100;
use Getopt::Long;
use Pod::Usage;

our $opfn="";
my $verbose;
my $ID2GSB ='./ID2GSB-101040.txt';
our $count2=0;

open (FILE,$ID2GSB)or die ("File ID2GSB-101040.txt not Found\n");

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
or die("Error in command line arguments\nUsage: perl transGSB2ID.pl [-o|outfile string] <input_file>\n");



if ($opfn eq ""){
$opfn="GSV2ID.out";
print "Output file:$opfn.csv \n";
}else{
print "Output file:$opfn.csv \n";
}
open OUT, "> $opfn.csv" or die ("[-] Error: Can't open or creat $opfn.txt\n");
print OUT "\"BNAID\",\"GSBID\",\"foldchange\",\"log2FoldChange\",\"stat\",\"pvalue\",\"padj\"\n";



our %GSB2ID;
while(my $seq2 = <FILE>){
	my @list=split/\t/,$seq2;
	#$atid=substr($at[0],0,9);
	my $GSBid= $list[1];
	$GSBid=~s/\n//;
	#say $GSBid;
	  $GSB2ID{$GSBid}=$list[0];
	 
}
close FILE;
print "Finished loading ID2GSB-101040.txt!\n";



while(defined(our $seq = <>)){
   if ($seq=~m/GSB/g) {
   my @tmp=split/\,/,$seq;
	 my $GSBtmpid=$tmp[0];
	 $GSBtmpid=~ s/\"//g;
	 #print "$GSBtmpid\n";
	 my $BNAID=$GSB2ID{$GSBtmpid};
	 print OUT "\"$BNAID\",$seq";
	 $count2++;
  }else{
  next;
  }
	
}

print "DONE! $count2 lines.\n";