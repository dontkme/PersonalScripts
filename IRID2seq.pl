#!/bin/env perl
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
our $opfn="IRseqOut";
my $verbose;
our $upstreml=200;
our $downstreml=200;
our $annot;
our $intronLT=400;
# our $sortID="ID";
# our $feature="mRNA";

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"u=i"=>\$upstreml,"d=i"=>\$downstreml, "i=s"=>\$annot,"s=i"=>\$intronLT)
or die("[-]Error in command line arguments
  Usage: perl IRID2seq.pl [options] <-i string|Intron ID file> <input Genome FASTA file>
    options:
    [-o string|outprefix Default: IRseqOut]
    [-s string|Intron Legnth cutoff default: 400]
    [-u int|upstream length Default: 200]
    [-d int|downstream length Default: 200]
	 
    Note: Get Sequences of IRID from genome v0.100 2022/03/20.\n");


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
open ANNOT, "< $annot" or die ("[-] Error: Can't open IRID file: $annot.\n");

open OUT, "> $opfn.in$intronLT.fa" or die ("[-] Error: Can't open or create $opfn.fa\n");
if ($upstreml >= 0){
    open OUTUP, "> $opfn.up$upstreml.fa" or die ("[-] Error: Can't open or create $opfn.up$upstreml.fa\n");
} else {
    die ("[-] Error: Upstream lenghth must be zero or greater\n");
}

if ($downstreml >= 0){
    open OUTDOWN, "> $opfn.down$downstreml.fa" or die ("[-] Error: Can't open or create $opfn.down$downstreml.fa\n");
} else {
    die ("[-] Error: Downstream lenghth must be zero or greater\n");
}

# open OUTALL, "> $opfn.u$upstreml.d$downstreml.fa" or die ("[-] Error: Can't open or creat $opfn.u$upstreml.d$downstreml.fa\n");

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

	if ($seq =~ m/^.*>/) {
	# $seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
    $seq=~ m/^\s*>\s*(\S+)\s*(.*)/;
	# print "$1\n";
	 $Chrname[$Chri]= $1;
	$Chri++;
	}else{
		$seq =~ s/\s//;
		$seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
	    #$Chrseq[$Chri-1] .=$seq;
     $Chrid2seq{$Chrname[$Chri-1]} .=$seq;
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
print "Running. Please wait for a minite.\n";
######### Finished loading.


#####################################
#Start main 
#####################################
print "Start loading input IRID.\n";
print "Output intron length < $intronLT: $opfn.in$intronLT.fa\n";
print "Output Down $downstreml: $opfn.down$downstreml.fa\n";
print "Output Up $upstreml: $opfn.up$upstreml.fa\n";

our @tmp;
our $faheadid="test";
our $annotcount=0;

while(defined(our $inrow = <ANNOT>)){
        $inrow =~ s/\"//g;
    if ($inrow =~ m/^\#/) {next;}
    if ($annotcount % 500 == 0){
        print "Dealed with $annotcount annotations.\n";
    }
    @tmp = split (/\//,$inrow);
    # say $inrow;
    # say scalar(@tmp);
    if (scalar(@tmp)==4){  
        $annotcount++;
        $inrow =~ s/\s//g;
        $faheadid = $inrow;
    }else{
        next;
    }
    # if ($tmp[2] ne $feature){next;}

    my @tmp2 =split (/:/,$tmp[3]);
    # say "tmp2:", scalar(@tmp2);
    my @tmp3 = split(/-/, $tmp2[1]);
    # say "tmp3:", scalar(@tmp3);
    # my $PM=$tmp2[2]; ## Plus Minus Strand.
    # $pm =~ s/\s//g; ## remove \n or \r\n.

    my $seqstartpos=$tmp3[0];
    my $seqendpos=$tmp3[1];
    my $plusminus=$tmp2[2];
    $plusminus =~ s/\s//g; ## remove \n or \r\n.
    my $seqchrid=$tmp2[0]; 
    # say $seqchrid, $seqstartpos, $seqendpos, $plusminus;

    my $ds=$seqendpos-$seqstartpos;
    if ($ds < $intronLT){

        if ($plusminus eq "+"){
            my $finalseq = substr($Chrid2seq{$seqchrid},$seqstartpos,$seqendpos-$seqstartpos);
            print OUT ">$faheadid"." IRL $ds\n";
            print OUT "$finalseq\n";
        
        }elsif($plusminus eq "-"){
            my $finalseq= TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos,$seqendpos-$seqstartpos));
            print OUT ">$faheadid"." IRL $ds\n";
            print OUT "$finalseq\n";
        }

    }else{

        if ($plusminus eq "+"){

            my $upfinalseq="";

            print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:".($seqendpos-$upstreml)."-$seqendpos $plusminus\n";
            $upfinalseq = substr($Chrid2seq{$seqchrid},$seqendpos-$upstreml,$upstreml); # seqendpos
            # print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:$seqstartpos..$seqendpos\n";
            print OUTUP "$upfinalseq\n";

            ### Down
            my $downfinalseq="";
            print OUTDOWN ">$faheadid.down$downstreml"."  $seqchrid:$seqstartpos-".($seqstartpos+$downstreml)." $plusminus\n";
            $downfinalseq = substr($Chrid2seq{$seqchrid},$seqstartpos,$downstreml);
            print OUTDOWN "$downfinalseq\n";


        }elsif($plusminus eq "-"){

            my $upfinalseq="";

            $upfinalseq = TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos,$upstreml));
            print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:$seqstartpos-".($seqstartpos+$upstreml)." $plusminus\n";
            print OUTUP "$upfinalseq\n";

            ### Down
            my $downfinalseq="";
            $downfinalseq = TRseq(substr($Chrid2seq{$seqchrid},$seqendpos-$downstreml,$downstreml));
            print OUTDOWN ">$faheadid.down$downstreml $seqchrid:", $seqendpos-$downstreml,'-' ,"$seqendpos $plusminus\n";
            print OUTDOWN "$downfinalseq\n";

        
        }

    }

  

}

close ANNOT;
print "All done. Dealed with $annotcount IRID.\n";
close OUT;
close OUTDOWN;
close OUTUP;

