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
our $opfn="getASOOut";
our $ASOprefix="";
my $verbose;
our $window=20;
our $step=6;
# our $annot;
# our $sortID="ID";
# our $feature="mRNA";

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"n=s" => \$ASOprefix,"verbose"=>\$verbose,"w=i"=>\$window,"s=i"=>\$step)
or die('[-]Error in command line arguments
  Usage: perl GetASOfromFA.pl [options] <input FASTA file>
    options:
    [-o string|outprefix Default: getASOOut]
    [-n string|ASO prefix: Default: "", use fasta header when NULL]
    [-w int|window length Default: 20]
    [-s int|step length Default: 6]
	 
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

open OUT, "> $opfn.w$window.s$step.fa" or die ("[-] Error: Can't open or create $opfn.w$window.s$step.fa\n");
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

open OUTALL, "> $opfn.w$window.s$step.txt" or die ("[-] Error: Can't open or create $opfn.w$window.s$step.txt\n");


####################Output files End###########

################
# Loading Genome.
################
if (not $ARGV[0]) {
	die ("[-] Error: Not find a input Genome FASTA file.\n");
}
our $loadingstarttime=time();
# print @ARGV;
# if (@ARGV eq "") {
#     die ("[-] Error: Not find a input genome FASTA file.\n");
#     }
print "Start loading genomeic sequences.\n";

our $Chri=0;
our @Chrname=();
our @Chrseq=();
our %Chrid2seq;
#@ARGV = qw#''  Not_Find_a_File#;
#say @ARGV;
#say $0;
while(defined(our $seq = <>)){

	if ($seq =~ m/^.*>/){
	# $seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
    $seq=~ m/^\s*>\s*(\S+)\s*(.*)/;
	print "$1\n";
	 $Chrname[$Chri]= $1;
	$Chri++;
	}else{
		$seq =~ s/\s//;
		$seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
	    $Chrseq[$Chri-1] .=$seq;
    #  $Chrid2seq{$Chrname[$Chri-1]} .=$seq;
	}
}

    # for (our $i=0;$i<$Chri-1;$i++){
    #     $Chrid2seq{$Chrname[$i]}=$Chrseq[$i];
    # }

#close SEQFILENAME;
our $loadingendtime=time();
print "$Chri Sequences\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;



our $starttime=time();
print "Running. Please wait for a minute.\n";
#####################################
#Start main 
#####################################
print "Start getting ASO.\n";
print "Window: $window\n";
print "Step: $step\n";
print "ASO prefix: $ASOprefix\n";
# print "ID: $sortID\n";
# print "Feature: $feature\n";
our $annotcount=0;
our $targetseqLen=0;
for (my $ni=0;$ni<$Chri;$ni++){
    our $seqall = $Chrseq[$ni];
    # say $seqall;
    my $ChrID=$Chrname[$ni];

    my $seqname_prefix="";
    my $seqname_NO="";
    if ($ASOprefix eq ""){
        $seqname_prefix=$ChrID;
    }else{
        my $trueTagetNO=$ni+1;
        $seqname_prefix="$ASOprefix.$trueTagetNO";
    }
    # $targetseqLen=0;
    $targetseqLen=length($seqall);

    my $finalASOseq="";
    my $lastend=$targetseqLen - $window;
    my $ASONO=0;

    for (my $start=0; $start<$lastend; $start=$start+$step){
        $ASONO++;
        $finalASOseq=substr($seqall, $start, $window);
        $finalASOseq=TRseq($finalASOseq);

        my $finalASOname="$seqname_prefix.$ASONO";
        my $Notes = "$ChrID.$start.$window.$step";

        print OUT ">$finalASOname\n$finalASOseq\n";
        print OUTALL "$finalASOname\t$finalASOseq\t$Notes\n";
        # print "$finalASOseq\n";
    }


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
