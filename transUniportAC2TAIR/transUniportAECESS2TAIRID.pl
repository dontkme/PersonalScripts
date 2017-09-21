#!/usr/bin/perl
use strict;
#use warnings;
use 5.0100;

my $tair10dfilename ='./Tair10.33.314352uniportcombine.txt';
open (FILE,$tair10dfilename);
if(!open FILE,"< $tair10dfilename"){
 die "File not Found\n";
}
my $count2=0;
our $atid="";
our %unip2tair;
while(my $seq2 = <FILE>){
	my @at=split/\t/,$seq2;
	#print "@at[1]\n";
	if (@at[1] eq "-"){
	    next;
	  }else{
	    my @unipac=split/;\s/,@at[1];
 	#$atid=substr($at[0],0,9);
	    $atid=$at[0];
	    #say $atid;
	  foreach my $i (@unipac){
	  #say $i;
	    $i=~s/\s//g;
	     #say $i;
	  
	    $unip2tair{$i}=$atid;
	    #say $unip2tair{$i};
	    }
	 
	}

}
close FILE;
say "Finished loading!";


my $dir =<STDIN>;
chomp $dir;
my $difffilename ="./inputuniport.txt";
open (DIFFFILENAME,$difffilename);
if(!open DIFFFILENAME,"< $difffilename"){
 die "File not Found\n";
}
my $count3=0;
open NEW,"> ./OUTuniportacees2TairID.txt";
#print NEW "test_id\tgene_id\ttype\tShort_description\tCurator_summary\tComputational_description\n";
while(my $diff = <DIFFFILENAME>){
	$count3++;
	my @ac=split/\|/,$diff;
	my $uniac=@ac[1];
	#say $uniac;
	our $tairid=$unip2tair{$uniac};
	#say $tairid;
	#say $count3;

	print NEW "$uniac\t$tairid\n"; 
	

	
}
close NEW;
say "DONE! $count3 ";