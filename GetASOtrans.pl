#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2023
# Get ASO Sequences from target fasta v1.0000 2023/11/03
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
our $opfn="getASOtransOut";
our $ASOprefix="";
my $verbose;
# our $window=20;
# our $step=6;
# our $annot;
# our $sortID="ID";
# our $feature="mRNA";

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"n=s" => \$ASOprefix,"verbose"=>\$verbose)
or die('[-]Error in command line arguments
  Usage: perl GetASOtrans.pl [options] <input FASTA file>
    options:
    [-o string|outprefix Default: getASOtransOut]
    [-n string|ASO prefix: Default: "", use fasta header when NULL]
	 
    Note: Get ASO Sequences from target fasta v1.0000 2023/11/03\n');

###################sub TRseq##########
 


sub TRseq($)
{
	my ($pinseq) = @_;
	#say $pinseq;
	my $pinseqtr = reverse $pinseq;
	#say $pinseqtr;
	 $pinseqtr =~ tr/ACGTacgt/TGCAtgca/;
	 #say $pinseqtr;
	 return  $pinseqtr;
}
##################TRseq End#############

####################Output files###########
# open ANNOT, "< $annot" or die ("[-] Error: Can't open annot file: $annot.\n");

open OUT, "> $opfn.transASO.fa" or die ("[-] Error: Can't open or create $opfn.transASO.fa\n");
open OUTRAW, "> $opfn.ASOraw.fa" or die ("[-] Error: Can't open or create $opfn.ASOraw.fa\n");
# if ($upstreml >= 0){
#     open OUTUP, "> $opfn.up$upstreml.fa" or die ("[-] Error: Can't open or create $opfn.up$upstreml.fa\n");
# } else {
#     die ("[-] Error: Upstream lenghth must be zero or greater\n");
# }

# if ($downstreml >= 0){
#     open OUTDOWN, "> $opfn.down$downstreml.fa" or die ("[-] Error: Can't open or create $opfn.down$downstreml.fa\n");
# } else {
#     die ("[-] Error: Downstream lenghth must be zero or greater\n");
# }

open OUTALL, "> $opfn.transASO.txt" or die ("[-] Error: Can't open or create $opfn.transASO.txt\n");


####################Output files End###########

################
# Loading Genome.
################
if (not $ARGV[0]) {
	die ("[-] Error: Not find a input ASO list file.\n");
}

# print @ARGV;
# if (@ARGV eq "") {
#     die ("[-] Error: Not find a input genome FASTA file.\n");
#     }
# print "Start loading genomeic sequences.\n";




our $starttime=time();
print "Running. Please wait for a minute.\n";
#####################################
#Start main 
#####################################
print "Start getting ASO.\n";
# print "Window: $window\n";
# print "Step: $step\n";
print "ASO prefix: $ASOprefix\n";
# print "ID: $sortID\n";
# print "Feature: $feature\n";
    our @tmp;
    our $faheadid="test";
    our $annotcount=0;
while(defined(our $inrow = <>)){

    if ($inrow =~ m/^\#/) {next;}

    if ($annotcount % 1000 == 0){
        print "Dealing $annotcount annotations.\n";
    }
    
    # chomp $inrow;
    $inrow =~ s/\R//g;
    # say $inrow;

    @tmp = split (/\t/,$inrow);
    my $ASOname=$tmp[0];
    my $ASOseq=$tmp[1];
    $ASOseq = uc $ASOseq;
    my $ASOnote=$tmp[2];
    my $ASOlength=length($ASOseq);
    my $ASOlast=$ASOlength-1;
    my $ASO2ndlast=$ASOlength-2;

    my @ASOnt= split(//, $ASOseq);

    $ASOnt[0]=  "/52MOEr"."$ASOnt[0]";
    my $ASOtransSeq=join("/*/i2MOEr", @ASOnt[0..$ASO2ndlast]);
    # say $ASOtransSeq;
    $ASOnt[$ASOlast]=  "$ASOnt[$ASOlast]"."/";
    $ASOtransSeq=join("/*/32MOEr",$ASOtransSeq,$ASOnt[$ASOlast]);
    # say $ASOtransSeq;

    my $finalline=join("\t",$ASOname,$ASOseq, $ASOtransSeq, $ASOnote);

    print OUTALL "$finalline\n";
    print OUT ">$ASOname $ASOnote\n$ASOtransSeq\n";
    print OUTRAW ">$ASOname $ASOnote\n$ASOseq\n";



    $annotcount++;
}   
  
  


print "All done. Dealed with $annotcount target sequence(s).\n";
close OUT;
close OUTALL;


######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;
