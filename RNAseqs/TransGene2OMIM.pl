#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2021
# Translate Gene Symbol to OMIM annotation v1.00 2021/10/05.
# hukaining@gmail.com

use strict;
# use warnings;
use 5.0100;

use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';

our $opfn="TransGene2OMIM.out";
my $verbose;
our $OMIM="";

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"i=s"=>\$OMIM)
or die("[-] Error in command line arguments
  Usage: perl TransGene2OMIM.pl [options] -i OMIM_genemap2_file <input gene name list file>
    options:
	 [-o string|output prefix Default: TransGene2OMIM.out]

    Note: Translate Gene Symbol to OMIM annotation v1.00 2021/10/05.\n");

### Loading OMIM annotations.
if (not @ARGV) {
	die ("[-] Error: Not find a gene list file as input.\n");
}

open ANNOT, "< $OMIM" or die ("[-] Error: Can't open OMIM_genemap2_file file: $OMIM.\n");
open OUT, "> $opfn.OMIMannot.txt" or die ("[-] Error: Can't open or create $opfn.OMIMannot.txt\n"); # Check output file writable.

our $loadingstarttime=time();

print "Start loading OMIM genemap2 file: $OMIM\n";
our $count1 = 0; # Count annotaion line numbers.
our $count2 = 0; # Count key Gene Symbol numbers.
our (%key2Chr, %key2StartPos, %key2EndPos, %key2CytoLoc, %key2ComputedCytoLoc, %key2MIMNum, %key2GeneSymbols,%key2GeneName,%key2ApprovedGeneSymbol,%key2EntrezID,%key2EnsemblID,%key2Comments,%key2Phenotypes,%key2MouseGeneSymbol);

while(defined(our $annotline=<ANNOT>)){
    if ($annotline =~ m/^#/){ # Skip '#' Lines.
        next;
    }

    $count1++;
    my @tmp=();
    @tmp = split("\t", $annotline);
    my @genesymbols = split(", ", $tmp[6]); # Get single Gene symbol.
    foreach my $skey (@genesymbols){
        $count2++;
        # say $skey;
        $key2Chr{$skey} = $tmp[0];
        $key2StartPos{$skey} = $tmp[1];
        $key2EndPos{$skey} = $tmp[2];
        $key2CytoLoc{$skey} = $tmp[3];
        $key2ComputedCytoLoc{$skey} = $tmp[4];
        $key2MIMNum{$skey} = $tmp[5];
        $key2GeneSymbols{$skey} = $tmp[6];
        $key2GeneName{$skey} = $tmp[7];
        $key2ApprovedGeneSymbol{$skey} = $tmp[8];
        $key2EntrezID{$skey} = $tmp[9];
        $key2EnsemblID{$skey} = $tmp[10];
        $key2Comments{$skey} = $tmp[11];
        $key2Phenotypes{$skey} = $tmp[12];
        chomp($tmp[13]);
        $key2MouseGeneSymbol{$skey} = $tmp[13];

    }

    if ($tmp[8] ne ""){
        my $skey = $tmp[8];
        if(!defined($key2MIMNum{$skey})){

            $count2++;
            # say $skey;
            $key2Chr{$skey} = $tmp[0];
            $key2StartPos{$skey} = $tmp[1];
            $key2EndPos{$skey} = $tmp[2];
            $key2CytoLoc{$skey} = $tmp[3];
            $key2ComputedCytoLoc{$skey} = $tmp[4];
            $key2MIMNum{$skey} = $tmp[5];
            $key2GeneSymbols{$skey} = $tmp[6];
            $key2GeneName{$skey} = $tmp[7];
            $key2ApprovedGeneSymbol{$skey} = $tmp[8];
            $key2EntrezID{$skey} = $tmp[9];
            $key2EnsemblID{$skey} = $tmp[10];
            $key2Comments{$skey} = $tmp[11];
            $key2Phenotypes{$skey} = $tmp[12];
            chomp($tmp[13]);
            $key2MouseGeneSymbol{$skey} = $tmp[13];

        }


    }


}

our $loadingendtime=time();
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;
print "$count1 annotations, $count2 gene symbols.\n";
print "====================================================================================\n";
# print "$key2GeneName{ADD1}\n";
# print "$key2GeneName{CRTC2}\n";

###### Start annotation.
$loadingstarttime=time();
print "Start annotaion.\n";
my $header = join("\t", "Chromosome\tGenomic Position Start\tGenomic Position End\tCyto Location\tComputed Cyto Location\tMIM Number\tGene Symbols\tGene Name\tApproved Gene Symbol\tEntrez Gene ID\tEnsembl Gene ID\tComments\tPhenotypes\tMouse Gene Symbol/ID");
# print OUT "$header\n";
our $count3=0;
our $count4=0;
while (defined(my $line = <>)){
    chomp($line);
    if ($line =~ m/^#/){ # Skip '#' Lines.
        next;
    }
    $count3++;
    my @input = split("\t",$line);

    if ($count3 == 1){  ### Print header. if input have multiple columns.

        my $columnumber = (scalar @input) - 1;
        print OUT "Gene symbol";
        for (my $i =0;$i<=$columnumber;$i++){
            print OUT "\t";
        }
        print OUT "$header\n";

    }

    my $tmpkey = uc $input[0]; # uppercase
    my $rest = join("\t",$key2Chr{$tmpkey}, $key2StartPos{$tmpkey}, $key2EndPos{$tmpkey}, $key2CytoLoc{$tmpkey}, $key2ComputedCytoLoc{$tmpkey}, $key2MIMNum{$tmpkey}, $key2GeneSymbols{$tmpkey},$key2GeneName{$tmpkey},$key2ApprovedGeneSymbol{$tmpkey},$key2EntrezID{$tmpkey},$key2EnsemblID{$tmpkey},$key2Comments{$tmpkey},$key2Phenotypes{$tmpkey},$key2MouseGeneSymbol{$tmpkey});
    if ($key2MIMNum{$tmpkey} ne ""){
        $count4++;
    }
    print OUT "$line\t$rest\n"; # print output.
}
$loadingendtime=time();
print "The annotation is complete!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;
print "$count3 input gene symbols, $count4 have OMIM annotations.\n";
print "Output file: $opfn.OMIMannot.txt\n";
