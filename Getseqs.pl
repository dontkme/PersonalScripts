#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2018
# Get Sequences from GFF with genome v1.0000 2018/09/09
# hukaining@gmail.com

use strict;
use warnings;
#use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="getseqsOut";
my $verbose;
our $upstreml=5000;
our $downstreml=5000;
our $annot;
our $sortID="ID";
our $feature="mRNA";

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"u=i"=>\$upstreml,"d=i"=>\$downstreml, "g=s"=>\$annot,"f=s"=>\$feature,"s=s"=>\$sortID)
or die("[-]Error in command line arguments
  Usage: perl EAHelitron [options] <-g string|GFF annoation file> <input FASTA file>
    options:
    [-o string|outprefix Default: getseqsOut]
    [-s string|Specify attribute type in GFF annotation for sorting. default: ID]
    [-f string|Specify feature type in GFF annotation.default: mRNA]
    [-u int|upstream length Default: 5000]
    [-d int|downstream length Default: 5000]
	 
    Note: Get Sequences from GFF with genome v1.0000 2018/09/09.\n");

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
open ANNOT, "< $annot" or die ("[-] Error: Can't open annot file: $annot.\n");

open OUT, "> $opfn.fa" or die ("[-] Error: Can't open or creat $opfn.fa\n");
if ($upstreml >= 0){
    open OUTUP, "> $opfn.up$upstreml.fa" or die ("[-] Error: Can't open or creat $opfn.up$upstreml.fa\n");
} else {
    die ("[-] Error: Upstream lenghth must be zero or greater\n");
}

if ($downstreml >= 0){
    open OUTDOWN, "> $opfn.down$downstreml.fa" or die ("[-] Error: Can't open or creat $opfn.down$downstreml.fa\n");
} else {
    die ("[-] Error: Downstream lenghth must be zero or greater\n");
}
####################Output files End###########

################
# Loading Genome.
################
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
	$seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
	print "$1\n";
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
#####################################
#Start main 
#####################################
print "Start loading input GFF.\n";
print "ID: $sortID\n";
print "Feature: $feature\n";


    

    our @tmp;
    our $faheadid="test";
    our $annotcount=0;
while(defined(our $inrow = <ANNOT>)){

    if ($inrow =~ m/^\#/) {next;}
    if ($annotcount % 1000 == 0){
        print "Dealed with $annotcount annotations.\n";
    }
    @tmp = split (/\t/,$inrow);
    # say $inrow;
    
    if ($tmp[2] ne $feature){next;}
    my @tmp2 =split (/\;/,$tmp[8]);

    foreach my $tmp3 (@tmp2){
        # if ($tmp3 =~ m/$sortID\=(\S+)/i){
        if ($tmp3 =~ m/$sortID\=(\S+)/i){
            $faheadid=$1;
            $annotcount++;
        }
        # } else {next;}
    }

    my $seqstartpos=$tmp[3];
    my $seqendpos=$tmp[4];
    my $plusminus=$tmp[6];
    my $seqchrid=$tmp[0];

    if ($plusminus eq "+"){
        my $finalseq="";
        $finalseq= substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1);
        print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
        print OUT "$finalseq\n";

        my $upfinalseq="";
        if ($seqstartpos<$upstreml){
            print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:1..$seqstartpos $plusminus\n";
            $upfinalseq = substr($Chrid2seq{$seqchrid},1,$seqstartpos);
        }else{
            print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:".($seqstartpos-$upstreml)."..$seqstartpos $plusminus\n";
            $upfinalseq = substr($Chrid2seq{$seqchrid},$seqstartpos-$upstreml,$upstreml);

        }
        # print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:$seqstartpos..$seqendpos\n";
        print OUTUP "$upfinalseq\n";

        my $downfinalseq="";
        $downfinalseq = substr($Chrid2seq{$seqchrid},$seqendpos,$downstreml);
        print OUTDOWN ">$faheadid.down$downstreml"."  $seqchrid:$seqendpos..".($seqstartpos+$downstreml)." $plusminus\n";
        print OUTDOWN "$downfinalseq\n";

    
    } elsif($plusminus eq "-"){
        my $finalseq= TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1));
        print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
        print OUT "$finalseq\n";

        my $downfinalseq="";
        if ($seqstartpos<$downstreml){
            print OUTDOWN ">$faheadid.down$downstreml"."  $seqchrid:1..$seqstartpos $plusminus\n";
            $downfinalseq = TRseq(substr($Chrid2seq{$seqchrid},0,$seqstartpos-1));
            print OUTDOWN "$downfinalseq\n";
        }else{
            print OUTDOWN ">$faheadid.down$downstreml"."  $seqchrid:".($seqstartpos-$downstreml)."..$seqstartpos $plusminus\n";
            $downfinalseq = TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos-$downstreml-1,$downstreml));
            print OUTDOWN "$downfinalseq\n";

        }
        # print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:$seqstartpos..$seqendpos\n";
        

        my $upfinalseq="";
        $upfinalseq = TRseq(substr($Chrid2seq{$seqchrid},$seqendpos,$upstreml));
        print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:$seqendpos..".($seqendpos+$upstreml)." $plusminus\n";
        print OUTUP "$upfinalseq\n";

    }


            
        
    

}

close ANNOT;
print "All done. Dealed with $annotcount annotations.\n";
close OUT;
close OUTDOWN;
close OUTUP;

######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;
