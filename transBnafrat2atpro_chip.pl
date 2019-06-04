#!/usr/bin/perl
use strict;
use warnings;
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
our (%at2type, %at2short, %at2sum, %at2detail);
my $count2=0;
while(my $seq2 = <FILE>){
	my @at=split/\t/,$seq2;
	#$atid=substr($at[0],0,9);
	our $atid=$at[0];
	  $at2type{$atid}=$at[1];
	 $at2short{$atid}=$at[2];
	  $at2sum{$atid}=$at[3];
	$_=$at[4];
	s/\n//;
	$at2detail{$atid}=$_;
}
close FILE;

my $bn2atidfile ='D:/transDB/Bnafrcds2tai10blastn1e-5.outfmt6';
open (FILE,$bn2atidfile) or die "File not Found\n";
our %bn2atid;
  $count2=0;
while(my $seq2 = <FILE>){
	my @at=split/\t/,$seq2;
	#$atid=substr($at[0],0,9);
	my $bnid=$at[0];
	  $bn2atid{$bnid}=$at[1];
	 #$at2short{$atid}=$at[2];
	  #$at2sum{$atid}=$at[3];
	#$_=$at[4];
	#s/\n//;
	#$at2detail{$atid}=$_;
}
close FILE;


say "Finished loading!";


my $dir =<STDIN>;
chomp $dir;
my $difffilename ="./ZK12PE.narrowPeak.overlap.bed.annot.txt";
open (DIFFFILENAME,$difffilename);
if(!open DIFFFILENAME,"< $difffilename"){
 die "File not Found\n";
}
my $count3=0;
open NEW,"> ./ZK12PE.narrowPeak.overlap.bed.annot.addTair10.txt";
print NEW "test_id\tgene_id\ttype\tShort_description\tCurator_summary\tComputational_description\n";
while(my $diff = <DIFFFILENAME>){
	$count3++;
	#say $count3;
	#say $diff;
	#while ($diff =~ /(Bna.+D)\t(AT.G[0-9]+\.[0-9])\t(.*)\n/igc){
	#while ($diff){

	  my @bn=split/\t/,$diff;
	#say $1;
	#say $n2a{$1};
	#my $nat=$n2nrf{$2};
	#say $2;
	$bn[21] =~s/\n//g;
	my $bnid=$bn[19];
	my $atid=$bn2atid{$bnid};
	my $orgseq= join("\t",@bn);
	print NEW "$orgseq\t$atid\t$at2type{$atid}\t$at2short{$atid}\t$at2sum{$atid}\t$at2detail{$atid}\n"; 
	

	#}
}
close NEW;
say "DONE!";
