#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2023
# Subset vcf from GTF regions v1.000 2023/11/17
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
# use Parallel::ForkManager;
our $MAX_processes=1;
# my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
use Array::Utils qw(:all);
use List::MoreUtils ':all';
# use Tie::Hash::Regex;

our $opfn="getseqsOut";
my $verbose;
# our $upstreml=5000;
# our $downstreml=5000;
# our $annot = "";
our $sortID = "gene_id";
our $feature = "exon";
our $inputVCF="";
our $dj = 50;

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
# GetOptions("o=s" => \$opfn, "p=i"=>\$MAX_processes, "d=i"=>\$dj, "g=s"=>\$annot,"f=s"=>\$feature,"s=s"=>\$sortID, "input|in=s"=>\$inputlist,"verbose"  => \$verbose)
GetOptions("o=s" => \$opfn, "d=i"=>\$dj,"f=s"=>\$feature,"s=s"=>\$sortID, "r|reference=s"=>\$inputVCF,"verbose"  => \$verbose)
or die("[-]Error in command line arguments
  Usage: perl GetTransIDsfromGTF [options]  -r <string|input indexed VCF file> <input GTF>
    options:
    [-o string|outprefix Default: getseqsOut]
    [-s string|Specify attribute type in GTF annotation for sorting. default: gene_id]
    [-f string|Specify feature type in GTF annotation.default: '']
       
	 
    Note: Subset vcf from GTF regions v1.0000 2023/11/17.\n");

# open ANNOT, "< $annot" or die ("[-] Error: Can't open annot file: $annot.\n");
# if (not @ARGV) {
# 	die ("[-] Error: Not find a input GTF file.\n");
# }


open OUT, "> $opfn.subset.txt" or die ("[-] Error: Can't open or create $opfn.subset.txt\n");

print "Start loading input GTF. (@ARGV)\n";
print "ID: $sortID\n";
print "Feature: $feature\n";
print "Indexed VCF: $inputVCF\n";

our @tmp=();
# our $faheadid="test";
# our $annotcount=0;
# our (%key2Chr,%key2type, %key2startpos, %key2endpos, %key2PM,%key2restwords,%key2line,%key2transid,%key2geneid,%key2exonid,%key2exonnumber,%SEexonid);
# our @codingtransid=();
# our (%trans_exonnumbers,%key2length, %trans_CDSstart, %trans_CDSend, %trans_PM);

our $transid ="";
our $geneid="";
our $EnsID="";
our $exonNumber="";
our $exonPM="";
our $exonStart="";
our $exonEnd="";
our $exonChr="";
our $tmpexonID="";
our $exonLen=0;
our $outputheaser="";

our $starttime=time();

while(defined (our $inrow=<>)){    ### Get transcripts id.
    
    if ($inrow =~ m/^#/){
        next;
    }

    $transid ="";
    $geneid="";
    $EnsID="";
    $exonNumber="";
    $exonPM="";
    $exonStart="";
    $exonEnd="";
    $exonChr="";
    $tmpexonID="";
    $exonLen=0;
    @tmp=();
    $outputheaser="";
    
    $inrow =~ s/\R//g;
    @tmp = split (/\t/,$inrow);

    $exonChr=$tmp[0];
    (my $decayChr = $exonChr)=~ s/chr//g;
    # say $decayChr;
    $exonStart=$tmp[3];
    $exonEnd=$tmp[4];
    $exonPM=$tmp[6];

    if($tmp[2] eq $feature){
    
        if ($inrow =~ m/transcript_id \"(\S+)\"/i){
           
            # push (@codingtransid, $1);
            $transid = $1;
        }

        if ($inrow =~ m/gene_id \"(\S+)\"/i){
            # push (@codingtransid, $1);
            $EnsID = $1;
        }

        if ($inrow =~ m/gene_name \"(\S+)\"/i){
            # push (@codingtransid, $1);
            $geneid = $1;
        }

        if ($inrow =~ m/exon_number (\S+);/i){
            # push (@codingtransid, $1);
            $exonNumber = $1;
        }

        $tmpexonID=join(":", $geneid, $exonChr,$exonStart, $exonEnd, $exonPM);

        $exonLen = $exonEnd - $exonStart+1;
        my $exonmod3 = $exonLen % 3;

        my $bcfTR = "$decayChr:$exonStart"."-". "$exonEnd";
        # say $bcfTR;
        $outputheaser=join("\t", $tmpexonID, $bcfTR, $geneid, $EnsID, $transid, $exonNumber, $exonLen, $exonmod3);
        # say $outputheaser;

        # my $bcfTres=system("bcftools view -H $inputVCF $bcfTR");# 2>/dev/null");
        my $bcfTres=`bcftools view -H -v snps $inputVCF $bcfTR`;# 2>/dev/null";
        # sleep(10);
        # print "$bcfTres\n";
        if ( $bcfTres ne 0){

            open TMPLINE ,"<", \$bcfTres;
            # print $bcfTres[0];
                my $linecount1=1;
            while (my $tmpbcfTresLine=<TMPLINE>){
                # say "line: $linecount1";
                $tmpbcfTresLine =~ s/\R//g;
                my $finalout=join("\t", $outputheaser, $tmpbcfTresLine);
                $finalout = "$finalout\n";
                print OUT $finalout;
                $linecount1++;
            }

            close TMPLINE;
        }



    } 
}

our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;