#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2021
# Compare CDD hitdata v1.00 2021/11/11
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
use Array::Utils qw(:all);
use List::MoreUtils ':all';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="cmpCDDhitdataOut";
our $n=0;
my $verbose;
# our $input;
# our $upstreml=5000;
# our $downstreml=5000;
# our $annot;
# our $sortID="ID";
# our $feature="mRNA";
GetOptions("o=s" => \$opfn,"n=i" => \$n,"verbose"=>\$verbose)
or die("[-]Error in command line arguments
  Usage: perl Getseqs [options] <-n Total input Q number> <input CDD hitdata file>
    options:
    [-o string|outprefix Default: getseqsOut]
    Note: Compare CDD hitdata v1.00 2021/11/11.\n");

#############

# open LIST, "< $input" or die ("[-] Error: Can't open input list file: $input.\n");
open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or create $opfn.txt\n");

####### Start loading the sequences.
if (not $ARGV[0]) {
	die ("[-] Error: Not find a input CDDhitdata file.\n");
}
our $loadingstarttime=time();
# print @ARGV;
# if (@ARGV eq "") {
#     die ("[-] Error: Not find a input genome FASTA file.\n");
#     }
print "Start loading genomeic sequences.\n";

# our $seqi=0;
# our @Chrname=();
# our @Chrseq=();
our %Name2seq;
our $seqname="";
our $seqnum=0;
our %Q;
our %Qname;
our $x=1;
our @QID;
for($x=1;$x<=$n;$x++){
    # my $Qarray="Q#$x";
    # say $Qarray;
    # $QID[$x]=$Qarray;
    @{$Q{$x}}=();
    # @{$Qarray}=();
    # @{$Qarray}=();
}

if($n<1){
    die ("[-] Error: Can't load Q# numbers\n");
}

while(defined(our $seq = <>)){
#while(our $seq = <SEQFILENAME>){
	if ($seq =~ m/^.*Q#/) {
        $seqnum++;
	    $seq=~ m/^(Q#(\d+)) - >(\S*)\t(.*)/; #1.52 new RE based on BioPerl
        # say $1;
        #   say $2;
        #   say $3;
        my @tmplines=split('\t',$seq);
        $Qname{$1}=$3;
        # $Q{$1}[];
        push(@{$Q{$2}},$tmplines[8]); ## CDD hitdata format[8];
        # push(@{$1},$tmplines[8]); ## CDD hitdata format[8];


        # my @seqnames=split('\|',$1);
        #  $seqname= $seqnames[1]; # 1 is Gencode pc_translations.fa headline pattern. 
        # #  say $seqname;
        # $seqnum++;
        # }else{
        # 	$seq =~ s/\s//;
        # 	# $seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
        #  $Name2seq{$seqname}.=$seq;
	}
}
our $Qnum=0;
$Qnum=keys %Q;
# $Qnum=$n;

if($Qnum>0){

    for(my $i=1;$i<$Qnum;$i+=2){
        my @a1=();
        my @a2=();
        if(@{$Q{$i}}){
             @a1=@{$Q{$i}};

        }else{
            @a1=("");
        }
        # print @a1;
        my $j=$i+1;
        if( @{$Q{$i+1}}){
            @a2=@{$Q{$j}};

        }else{
            @a2=("");
        }
        # print @a2;
        my @ma1=array_minus(@a1,@a2);
        my @ma2=array_minus(@a2,@a1);
        my $ma1all=join(",",@ma1);
        my $ma2all=join(",",@ma2);
        print OUT "$i\t$j\t$ma1all\t$ma2all\n";


    }

}else{
    die ("[-] Error: Can't load Q# lines\n");
}

close OUT;
our $loadingendtime=time();
print "$seqnum Sequences\n";
print "Q: $Qnum\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;

#### Loading input list of Transcript IDs. fomat 1st name\t 2ed name\t gene name\n
# print "Start loading input list.\n";
# our $inputlineNum=0;
# our $outputlineNum=0;
# while(defined(our $listline=<LIST>)){
# 	$inputlineNum++;
# 	$listline =~ s/\R//g;
# 	my @tmplistline=split("\t", $listline);
# 	my $firstName=$tmplistline[0];
# 	my $secondName=$tmplistline[1];
# 	my $genename=$tmplistline[2];

# 	if(defined($firstName) && defined($secondName)){
# 		$outputlineNum++;
# 		print OUT ">$firstName:$genename\n";
# 		print OUT "$Name2seq{$firstName}\n";
# 		print OUT ">$secondName:$genename\n";
# 		print OUT "$Name2seq{$secondName}\n";
		
# 	}
# }
# print "Input list lines: $inputlineNum\nOutput: $outputlineNum\n";
# close LIST;
close OUT;
