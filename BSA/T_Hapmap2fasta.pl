#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2019
# Trans T_HAPMAP to seprate fata fomat v1.0000 2019/09/16
# hukaining@gmail.com
#
use strict;
use warnings;
#use 5.0100;
use Getopt::Long;

use re 'eval';
our $opfn="Allsnp.fa";
our $sep="\t";
my $verbose;

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
or die("Error in command line arguments\nUsage: perl T_Hapmap2fasta.pl [-o outfile_name] <input_T_hapmap>\n");


if ($opfn eq ""){
	$opfn="Allsnp.fa";
	print "Output file:$opfn\n";
}else{
	print "Output file:$opfn\n";
}
open OUT, "> $opfn" or die ("[-] Error: Can't open or creat $opfn\n");


my $count=0;

our $gene="gene";
while(defined(our $seq = <>)){
	#chomp $seq;
	#print $seq;
	#@snp=sep($seq,"\t");
	my @snp=split/\t/,$seq;
	my @snp2=splice(@snp,1);
	#while ($seq =~ /^(.+)\t(GO\:\d{7})(\s)?.*$/g){
	#say $gene;
	#say $1;
	#say $2;
	
	
		print OUT ">$snp[0]\n";
		print OUT @snp2;
		  #$gene = "$1";
		  $count++;
		  #say $gene; 
	
	
		
	#}	
	
	}

close OUT;
print "DONE! $count samples\n";