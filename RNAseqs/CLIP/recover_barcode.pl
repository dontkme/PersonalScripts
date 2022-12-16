#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
use re 'eval';
my $verbose;

our $headbarcode_len = 5;
our $opfn = "Out.recover.fq";

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"h=i"=>\$headbarcode_len)
or die("[-]Error in command line arguments
  Usage: perl recover_barcode.pl [options] <input fastq file> 
    or
  zcat <input fastq.gz file> |perl recover_barcode.pl | gzip > <output file name>
    options:
     [-o string|output Default: Out.recover.fq]
     [-h int|Length of head barcode in 1st line Default: 5]

    Note: Recover 1st line barcode to sequence and quality line (for CLIP CTK eCLIP pipeline) v1.00 2022/12/16.\n");


# print "Barcode length: $headbarcode_len\n";
################
# Start loading
################

our $loadingstarttime=time();
# print "\nStart loading raw fastq.\n";

our $countall=0;
our $line_n=0;
our $tmpbarcode="";
our @tmpheadline=();
our $tmpbarcodeLen=0;
our $finalheadline="";
our $final2edline="";
our $tmpQ="IIIII";
our $finalQ="";

while(defined(our $seq = <>)){
    $seq =~ s/\s$//;
    $line_n++;

    if ($line_n % 4 == 1){
        $countall++;
        @tmpheadline = split(/:/, $seq);

        $tmpbarcode=$tmpheadline[0];
        $tmpbarcode =~ s/@//;
        $tmpbarcodeLen = length($tmpbarcode);

        # print "$tmpbarcode: $tmpbarcodeLen\n";
        
        $tmpheadline[1]="@"."$tmpheadline[1]";
        $finalheadline = join(":", @tmpheadline[1 .. $#tmpheadline]);
        print "$finalheadline\n";
    }


    if ($line_n % 4 == 2){
       
        $final2edline = "$tmpbarcode"."$seq";
        print "$final2edline\n";

    }

    if ($line_n % 4 == 3){
       
        print "$seq\n";

    }

    if ($line_n % 4 == 0){
       
        $tmpQ= "I" x $tmpbarcodeLen;
        $finalQ = "$tmpQ"."$seq";
        print "$finalQ\n";

    }






}





our $loadingendtime=time();
# print "$countall Sequences\n";
# print "Finished loading!\n";
# printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;