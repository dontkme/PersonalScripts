#!/usr/bin/perl
#use strict;
#use warnings;
use 5.0100;
#my $filename ='./merged_asm/merged.gtf';
#open (FILENAME,$filename);
#if(!open FILENAME,"< $filename"){
# die "File not Found\n";
#}
#my $count1=0;
#while(my $seq = <FILENAME>){
#	$count1++;
#	#say $count1;
#	#say $seq;
#	while ($seq =~ /(XLOC_.+)\"; transcript_id(.*)oId.+\"(.+)\";.+nearest_ref \"(.+)\"; class/igc){
#	#say $1;
#	#say $4;
#	 $n2oid{$1} = $3;
#	 $n2nrf{$1} = $4;
#
#	}
#}
## close FILENAME;

my $tair10dfilename ='D:/transDB/TAIR10_functional_descriptions_20160930.txt';
open (FILE,$tair10dfilename);
if(!open FILE,"< $tair10dfilename"){
 die "File not Found\n";
}
my $count2=0;
while(my $seq2 = <FILE>){
	my @at=split/\t/,$seq2;
	#$atid=substr($at[0],0,9);
	$atid=$at[0];
	  $at2type{$atid}=$at[1];
	 $at2short{$atid}=$at[2];
	  $at2sum{$atid}=$at[3];
	$_=$at[4];
	s/\n//;
	$at2detail{$atid}=$_;
}
close FILE;
say "Finished loading!";


my $dir =<STDIN>;
chomp $dir;
my $difffilename ="./Bnafrcds2tai10blastn1e-5.outfmt6";
open (DIFFFILENAME,$difffilename);
if(!open DIFFFILENAME,"< $difffilename"){
 die "File not Found\n";
}
my $count3=0;
open NEW,">> ./Bnafrcds2tai10blastn1e-5.detail.txt";
print NEW "test_id\tgene_id\ttype\tShort_description\tCurator_summary\tComputational_description\n";
while(my $diff = <DIFFFILENAME>){
	$count3++;
	#say $count3;
	#say $diff;
	while ($diff =~ /(Bna.+D)\t(AT.G[0-9]+\.[0-9])\t(.*)\n/igc){
	#say $1;
	#say $n2a{$1};
	#my $nat=$n2nrf{$2};
	#say $2;
	print NEW "$1\t$2\t$at2type{$2}\t$at2short{$2}\t$at2sum{$2}\t$at2detail{$2}\n"; 
	

	}
}
close NEW;
say "DONE!";
