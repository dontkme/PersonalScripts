#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2021
# Get NMD transid from GTF with genome v1.0000 2021/03/29
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
use Array::Utils qw(:all);
use List::MoreUtils ':all';
# use Tie::Hash::Regex;

our $opfn="getseqsOut";
my $verbose;
# our $upstreml=5000;
# our $downstreml=5000;
our $annot = "";
our $sortID = "gene_id";
our $feature = "";
our $inputlist;
our $dj = 50;

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn, "d=i"=>\$dj, "g=s"=>\$annot,"f=s"=>\$feature,"s=s"=>\$sortID, "input|in=s"=>\$inputlist,"verbose"  => \$verbose)
or die("[-]Error in command line arguments
  Usage: perl GetTransIDsfromGTF [options] -g <string|GTF annoation file> -in <string|input rMATs type list> <input FASTA file(s)>
    options:
    [-o string|outprefix Default: getseqsOut]
    [-s string|Specify attribute type in GFF annotation for sorting. default: gene_id]
    [-f string|Specify feature type in GFF annotation.default: '']
       
	 
    Note: Get Transid from GTF with genome v1.0000 2021/03/29.\n");

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
########################################




####################Output files###########
open ANNOT, "< $annot" or die ("[-] Error: Can't open annot file: $annot.\n");

open OUT, "> $opfn.annotseq.fa" or die ("[-] Error: Can't open or create $opfn.annotseq.fa\n");

open INPUTLIST, "< $inputlist" or die ("[-] Error: Can't open list file: $inputlist.\n");

open OUTSEQ, "> $opfn.outseq.txt" or die ("[-] Error: Can't open or create $opfn.outseq.txt\n");

# if ($upstreml >= 0){
#     open OUTUP, "> $opfn.up$upstreml.fa" or die ("[-] Error: Can't open or creat $opfn.up$upstreml.fa\n");
# } else {
#     die ("[-] Error: Upstream lenghth must be zero or greater\n");
# }

# if ($downstreml >= 0){
#     open OUTDOWN, "> $opfn.down$downstreml.fa" or die ("[-] Error: Can't open or creat $opfn.down$downstreml.fa\n");
# } else {
#     die ("[-] Error: Downstream lenghth must be zero or greater\n");
# }

# open OUTALL, "> $opfn.u$upstreml.d$downstreml.fa" or die ("[-] Error: Can't open or creat $opfn.u$upstreml.d$downstreml.fa\n");


####################Output files End###########

################
# Loading Genome.
################
if (not @ARGV) {
	die ("[-] Error: Not find a input Genome FASTA file.\n");
}
our $loadingstarttime=time();
# print @ARGV;
# if (@ARGV eq "") {
#     die ("[-] Error: Not find a input genome FASTA file.\n");
#     }
print "Start loading genomeic sequence file(s): @ARGV \n";

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
print "====================================================================================\n";


our $starttime=time();
print "Running. Please wait for a minite.\n";

################sub get sequence########

sub getseq
{
    my $sortkey= $_[0];
    my $finalseq="";
    my $seqchrid=$_[1];
    my $seqstartpos=$_[2];;
    my $seqendpos=$_[3];;

    if ($_[4] eq "+"){
            
            # my $finalseq="";
            $finalseq= substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1);
            # print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
            # print OUT ">$sortkey\n";
            # print OUT "$finalseq\n";

        }elsif($_[4] eq "-"){
               
                # my $finalseq="";
                $finalseq= TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1));
                # print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
                # print OUT ">$sortkey\n";
                # print OUT "$finalseq\n";

        }
        return $finalseq
}



#####################################
#Start main 
#####################################
print "Start loading input GTF. ($annot)\n";
print "ID: $sortID\n";
print "Feature: $feature\n";


    

    our @tmp;
    our $faheadid="test";
    our $annotcount=0;
    our (%key2Chr,%key2type, %key2startpos, %key2endpos, %key2PM,%key2restwords,%key2line,%key2transid,%key2geneid,%key2exonid,%key2exonnumber,%SEexonid);
    our @codingtransid=();
    our (%trans_exonnumbers,%key2length, %trans_CDSstart, %trans_CDSend, %trans_PM);



while(defined (our $inrow=<ANNOT>)){
    @tmp=();
    @tmp = split (/\t/,$inrow);

    if($tmp[2] eq "CDS"){
        # say $inrow;

        if ($inrow =~ m/transcript_id \"(\S+)\"/i){

            # say $inrow;
            # say $1;

            push (@codingtransid, $1);

        }

        } else {

        next;
    }


}

@codingtransid = uniq(@codingtransid); #use List::MoreUtils ':all'; uniq

# print "@codingtransid\n";

print scalar(@codingtransid)." protein-coding transcript id(s)\n";

our $grep_codingtransids = join '|', @codingtransid; 
# my @matched_results = grep { /$grep_pos/ } @exome; 
our %allcodingtransids = map {$_=>1} @codingtransid;

# sleep(2);

open ANNOT, "< $annot" or die ("[-] Error: Can't open annot file: $annot.\n");

    
while(defined(our $inrow = <ANNOT>)){

    @tmp=();

    if ($inrow =~ m/^\#/) {next;}
    if ($annotcount > 1 and $annotcount % 100000 == 0){
        print "Deal with $annotcount annotations.\n";
    }
    


            @tmp = split (/\t/,$inrow);
            # say $inrow;
            
            # if ($tmp[2] ne $feature){next;} # Could #
            my @tmp2 =split (/\;/,$tmp[8]);

            our ($gene_id, $transcript_id, $exon_number, $exon_id, $gene_name,$exon_numberf3)=("","",0,"","","000");

            foreach my $tmp3 (@tmp2){
                # if ($tmp3 =~ m/$sortID\=(\S+)/i){
                # if ($tmp3 =~ m/$sortID\=(\S+)/i){
                #     $faheadid=$1;
                #     $annotcount++;
                # }
                # } else {next;}
                if ($tmp3 =~ m/gene_id \"(\S+)\"/i){
                    $gene_id = $1;
                    }elsif($tmp3 =~ m/transcript_id \"(\S+)\"/i){
                        $transcript_id = $1;            
                    }elsif($tmp3 =~ m/exon_number \"(\d+)\"/i){
                        $exon_number = $1 + 0;
                        # $exon_numberf3 = printf "%3.0f",$exon_number;
                        $exon_numberf3 = sprintf "%03d",$exon_number;
                    }elsif($tmp3 =~ m/exon_id \"(\S+)\"/i){
                        $exon_id = $1;
                    }elsif($tmp3 =~ m/gene_name \"(\S+)\"/i){
                        $gene_name = $1;
                }         
                  
            }

    # foreach my $tmpcodingtransid (@codingtransid){
        # if (grep {$transcript_id} @codingtransid){ 
        # if ($transcript_id =~m/$grep_codingtransids/){ 
        # if ($grep_codingtransids =~ m/$transcript_id/){ 
        # if ($transcript_id ~~ @codingtransid){ 
        if (exists $allcodingtransids{$transcript_id}){ 
            # say $transcript_id;

            my $seqstartpos=$tmp[3];
            my $seqendpos=$tmp[4];
            my $plusminus=$tmp[6];
            my $seqchrid=$tmp[0];
            my $annottype=$tmp[2];

            our $key1 = "$transcript_id.$gene_id.$annottype.exon.$exon_numberf3";
            our $key = "$transcript_id.$gene_id.$annottype.ex.$exon_numberf3.$seqchrid:$seqstartpos..$seqendpos:$plusminus";
            $annotcount++;
            # print "$key\n";
            $key2Chr{$key}=$seqchrid;
            $key2type{$key}=$annottype;
            $key2startpos{$key}=$seqstartpos;
            $key2endpos{$key}=$seqendpos;
            $key2PM{$key}=$plusminus;
            $key2restwords{$key}=$tmp[8];
            $key2line{$key}=$inrow;
            $key2transid{$key} = $transcript_id;
            $key2geneid{$key} = $gene_id;
            $key2exonnumber{$key} = $exon_number;   
            $key2exonid{$key} = $exon_id;   

            $key2length{$key} = abs($seqendpos-$seqstartpos+1);
            $trans_PM{$transcript_id} = $plusminus;

            # if ($exon_number > $trans_exonnumbers{$transcript_id} or !defined($trans_exonnumbers{$transcript_id})){
            if (!defined($trans_exonnumbers{$transcript_id}) or $exon_number > $trans_exonnumbers{$transcript_id}){
                $trans_exonnumbers{$transcript_id} =$exon_number;
            }

            # if 
                    if ($annottype eq "CDS"){
                      
                        if (!defined($trans_CDSstart{$transcript_id} ) or ($exon_number < $trans_CDSstart{$transcript_id}) ){
                            $trans_CDSstart{$transcript_id} = $exon_number;
                        }

                       
                        if ( !defined($trans_CDSend{$transcript_id}) or ($exon_number > $trans_CDSend{$transcript_id})){
                            $trans_CDSend{$transcript_id} = $exon_number;
                        }
                    }

        }else{
            next;
        }
    # }


}
print "Finished loading GTF. Loaded $annotcount annotations.\n";


close ANNOT;
print "All done. Dealed with $annotcount annotations.\n";

# while ((my $key, my $value) = each %key2PM){  
# print "$key => $value\n";  
# }  

# foreach our $sortkey (sort keys %key2PM) {  
#  print "$sortkey => $key2PM{$sortkey}\n" ;  
# };  

our @sortkey = (sort keys %key2PM); # All-key array. @sortkey !!!
our @sortCDSkey = grep(/CDS/, @sortkey);
print scalar(@sortCDSkey)." CDS annotations.\n";# 2021/04/06
our @sortexonkey = grep(/exon/, @sortkey);
print scalar(@sortexonkey)." exon annotations.\n";


our (%CDSseqs, %exonseqs, %CDSstartPos, %CDSendPos, %exonstartPos, %exonendPos, %CDSsLength, %exonsLength, %accCDSlength, %reverseaccCDSlength);

foreach my $sortkey(@sortkey){
    # say $sortkey;
    # if ($sortkey =~m/exon/i){
    # if ($sortkey =~m/$feature/i){
    if ($sortkey =~ m/CDS|exon/){
        my $finalseq="";
        my $seqchrid=$key2Chr{$sortkey};
        my $seqstartpos=$key2startpos{$sortkey};
        my $seqendpos=$key2endpos{$sortkey};
        my $PM=$key2PM{$sortkey};
        $finalseq = getseq($sortkey,$seqchrid,$seqstartpos,$seqendpos,$PM);

        my $transid = $key2transid{$sortkey};  #make exon/CDS seq hash.
        my $transtype = $key2type{$sortkey};
        my $transexonno = $key2exonnumber{$sortkey};
        my $transexonstartpos = $key2startpos{$sortkey};
        my $transexonendpos = $key2endpos{$sortkey};
        my $transexonLength = $key2length{$sortkey};

        # my $exontypeseqhashname = "$transid.$transtype";
        # our $$exontype{"$transid"}[$transexonno] = $finalseq;
        if ($transtype eq "CDS"){
            $CDSseqs{$transid}[$transexonno] = $finalseq;
            $CDSstartPos{$transid}[$transexonno] = $transexonstartpos;
            $CDSendPos{$transid}[$transexonno] = $transexonendpos;
            $CDSsLength{$transid}[$transexonno] = $transexonLength;

        } elsif ($transtype eq "exon"){
            $exonseqs{$transid}[$transexonno] = $finalseq;
            $exonstartPos{$transid}[$transexonno] = $transexonstartpos;
            $exonendPos{$transid}[$transexonno] = $transexonendpos;
            $exonsLength{$transid}[$transexonno] = $transexonLength;
        } 


        print OUT ">$sortkey\n";
        print OUT "$finalseq\n";
    #  
    }
}
# say @sortkey;







close OUT;

######################## condon2aa ############################
#
my(%genetic_code) = (
   
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
);

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{

            # die "Bad codon '$codon'!!\n";
            return "*";
    }
}

sub seq2aa1st {

    my $dna= $_[0];
    my $protein1="";

    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
    my $codon = substr($dna,$i,3);
        $protein1 .= codon2aa($codon);
    
    }
     return $protein1;
}

sub seq2aa2ed {

    my $dna= $_[0];
    my $protein2="";

    for(my $i=1; $i < (length($dna) - 2) ; $i += 3) {
    my $codon = substr($dna,$i,3);
       $protein2 .= codon2aa($codon);
    
    }
     return $protein2;
}

sub seq2aa3rd {

    my $dna= $_[0];
    my $protein3="";

    for(my $i=2; $i < (length($dna) - 2) ; $i += 3) {
    my $codon = substr($dna,$i,3);
       $protein3 .= codon2aa($codon);
    
    }
     return $protein3;
}
#
####################### codon2aa ###############################



# open INPUTLIST, "< $inputlist" or die ("[-] Error: Can't open list file: $inputlist.\n"); # moved to head

# open OUTSEQ, "> $opfn.outseq.txt" or die ("[-] Error: Can't open or creat $opfn.outseq.txt\n");

while(defined(our $inputline = <INPUTLIST>)){
    if ($inputline =~ m/^\#/) {next;}
    # if ($inputline > 1 and $annotcount % 1000 == 0){
    #         print "Dealed with $annotcount annotations.\n";
    #     }
    $inputline =~ s/\r\n// ;
        @tmp = split (/\t/,$inputline);
        #  say $inputline;
        
        # if ($tmp[2] ne $feature){next;}
        my $genename = $tmp[0];
        my @tmp2 =split (/@/,$tmp[1]);
        my @tmpSE =split (/\:/,$tmp2[0]);
        my @tmpUS =split (/\:/,$tmp2[1]);
        my @tmpDS =split (/\:/,$tmp2[2]);
        my ($SEchr, $SEstartpos, $SEendpos, $SEPM) = split (/\:/,$tmp2[0]);
        # say $SEchr;
        # say $SEstartpos;
        # say $SEendpos;
        # say $SEPM;
        my ($USchr, $USstartpos, $USendpos, $USPM) = split (/\:/,$tmp2[1]);
        my ($DSchr, $DSstartpos, $DSendpos, $DSPM) = split (/\:/,$tmp2[2]);

        my $USseq = getseq("$genename.UE",$USchr, $USstartpos+1, $USendpos, $USPM);
        my $SEseq = getseq("$genename.SE",$SEchr, $SEstartpos+1, $SEendpos, $SEPM);
        my $DSseq = getseq("$genename.DE",$DSchr, $DSstartpos+1, $DSendpos, $DSPM);
        my $Oriseq = "$USseq"."$SEseq"."$DSseq";
        my $SEedseq = "$USseq"."$DSseq";
        my $Oriseq1stFrameAA = seq2aa1st($Oriseq);
        my $SEedseq1stFrameAA = seq2aa1st($SEedseq);
        my $Oriseq1stFrameAACount = ($Oriseq1stFrameAA =~ tr/_/_/);
        # say $Oriseq1stFrameAACount;
        my $SEedseq1stFrameAACount =($SEedseq1stFrameAA =~ tr/_/_/);
        # say $SEedseq1stFrameAACount;
        my $Oriseq2edFrameAA = seq2aa2ed($Oriseq);
        my $SEedseq2edFrameAA = seq2aa2ed($SEedseq);
        my $Oriseq2edFrameAACount= ($Oriseq2edFrameAA =~ tr/_/_/);
        my $SEedseq2edFrameAACount =($SEedseq2edFrameAA =~ tr/_/_/);
        my $Oriseq3rdFrameAA = seq2aa3rd($Oriseq);
        my $SEedseq3rdFrameAA = seq2aa3rd($SEedseq);
        my $Oriseq3rdFrameAACount= ($Oriseq3rdFrameAA =~ tr/_/_/);
        my $SEedseq3rdFrameAACount =($SEedseq3rdFrameAA =~ tr/_/_/);

        print OUTSEQ "$inputline\t$USseq\t$SEseq\t$DSseq\t$Oriseq\t$SEedseq\t$Oriseq1stFrameAA\t$SEedseq1stFrameAA\t$Oriseq1stFrameAACount\t$SEedseq1stFrameAACount\t$Oriseq2edFrameAA\t$SEedseq2edFrameAA\t$Oriseq2edFrameAACount\t$SEedseq2edFrameAACount\t$Oriseq3rdFrameAA\t$SEedseq3rdFrameAA\t$Oriseq3rdFrameAACount\t$SEedseq3rdFrameAACount\n";



}

close OUTSEQ;
close INPUTLIST;
# close OUTDOWN;
# close OUTUP;
# use Array::Utils qw(:all);
# use List::MoreUtils ':all';

open INPUTLIST, "< $inputlist" or die ("[-] Error: Can't open list file: $inputlist.\n");
open OUTTRANSID, "> $opfn.outtransid.txt" or die ("[-] Error: Can't open or create $opfn.outtransid.txt\n");
open OUTSETRANSID, "> $opfn.outSEtransid.txt" or die ("[-] Error: Can't open or create $opfn.outSEtransid.txt\n");
open OUTUSTRANSID, "> $opfn.outUStransid.txt" or die ("[-] Error: Can't open or create $opfn.outUStransid.txt\n");
open OUTDSTRANSID, "> $opfn.outDStransid.txt" or die ("[-] Error: Can't open or create $opfn.outDStransid.txt\n");
open OUTUSDSTRANSID, "> $opfn.outUSDStransid.txt" or die ("[-] Error: Can't open or create $opfn.outUSDStransid.txt\n");
open OUTUSSEDSTRANSID, "> $opfn.outUSSEDStransid.txt" or die ("[-] Error: Can't open or create $opfn.outUSSEDStransid.txt\n");


# tie %key2transid, 'Tie::Hash::Regex';


while(defined(our $inputline = <INPUTLIST>)){
    if ($inputline =~ m/^\#/) {next;}
    # if ($inputline > 1 and $annotcount % 1000 == 0){
    #         print "Dealed with $annotcount annotations.\n";
    #     }
    $inputline =~ s/\r\n// ;

        my @tmp = split (/\t/,$inputline);
    
        my $genename = $tmp[0];
        my @tmp2 =split (/@/,$tmp[1]);
        my @tmpSE =split (/\:/,$tmp2[0]);
        my @tmpUS =split (/\:/,$tmp2[1]);
        my @tmpDS =split (/\:/,$tmp2[2]);
        # my ($SEchr, $SEstartpos, $SEendpos, $SEPM) = split (/\:/,$tmp2[0]);
        my $SEstartpos = $tmpSE[1]+1;
        my $USstartpos = $tmpUS[1]+1;
        my $DSstartpos = $tmpDS[1]+1;
        if ($tmpSE[3] eq ""){next;} #pass Null line.
        my $SEpos="$tmpSE[0]\:$SEstartpos\[.\]\{2\}$tmpSE[2]\:\\$tmpSE[3]";
        my $USpos="$tmpUS[0]\:$USstartpos\\.\\.$tmpUS[2]\:\\$tmpUS[3]";
        my $DSpos="$tmpDS[0]\:$DSstartpos\\.\\.$tmpDS[2]\:\\$tmpDS[3]";

        # my $SEposre="qr/$tmpSE[0]\:$SEstartpos\\.\\.$tmpSE[2]\\:$tmpSE[3]/";
        # my $USposre="qr/$tmpUS[0]\:$USstartpos\\.\\.$tmpUS[2]:\\$tmpUS[3]/";
        # my $DSposre="qr/$tmpDS[0]\:$DSstartpos\\.\\.$tmpDS[2]:\\$tmpDS[3]/";

        #  say $SEpos;
        #  say $SEposre;
        #  say $USpos;
        # say $USposre;
        #  say $DSpos;
        # say $DSposre;

        our (@pos1, @pos2, @pos3, @pos1trans)= ((),(),()); 
        # our (@pos1, @pos2, @pos3)= (); 


        # my @USposlist   = tied(%key2transid)->FETCH($USpos);
        

        foreach my $keypos1(@sortkey){
            if ($keypos1 =~ m/$USpos/){
                # say $keypos1;
                 my $tmptransid = $key2transid{$keypos1};
                 
                 push(@pos1, $tmptransid);

            }
        }

        @pos1 = uniq(@pos1);
        # @pos1 = uniq(@USposlist);

        # print "pos1: @pos1\n";

        # foreach my $UPposkey(@sortkey){
        #     my $tmptransid = $key2transid{$UPposkey}; 
        # }

        # my @DSposlist   = tied(%key2transid)->FETCH($DSpos);

        foreach my $keypos2(@sortkey){
            if ($keypos2 =~ m/$DSpos/){
                my $tmptransid = $key2transid{$keypos2};
                 push(@pos2, $tmptransid);
            }
        }

        @pos2 = uniq(@pos2);
        # @pos2 = uniq(@DSposlist);
        # print "pos2: @pos2\n";    

        # say $SEpos;

        # my @SEposlist   = tied(%key2transid)->FETCH($SEpos);
        # print @SEposlist."\n";

        foreach my $keypos3(@sortkey){
            if ($keypos3 =~ m/$SEpos/){
                my $tmptransid = $key2transid{$keypos3};
                #   $SEexonid{$tmptransid} = $key2exonnumber{$keypos3};
                 push(@pos3, $tmptransid);
            }
        }



        @pos3 = uniq(@pos3);
        # my @pos4 = uniq(@SEposlist);
        # print "pos3: @pos3\n";
        # print "pos4: @pos4\n";
        

        my @USDStransid = intersect(@pos1,@pos2);
        # print "US_DS_intersect_transids: @USDStransid\n";
        my @noSEtransid = array_minus(@USDStransid,@pos3);
        # print "NoSEtransid_array_minus: @noSEtransid\n";
        my @USSEDStransid = intersect(@USDStransid,@pos3);
        my @SEnoUStransid = array_minus(@pos3,@pos1);
        my @onlySEtransid = array_minus(@SEnoUStransid,@pos2);
        my @uniqueSEtransid = @pos3;
        # print "Unique_SE_transid: @uniqueSEtransid\n";


        # say $SEchr;
        # say $SEstartpos;
        # say $SEendpos;
        # say $SEPM;
        print OUTTRANSID "$inputline\t@uniqueSEtransid\t@pos1\t@pos2\t@USDStransid\t@noSEtransid\n";
        
        # foreach my $tmptransid(@uniqueSEtransid){
        #     my $setransexonid = $SEexonid{$tmptransid};
        #     print OUTSETRANSID "$inputline\t$tmptransid\t$setransexonid\n"
        # }

        foreach my $tmptransid(@onlySEtransid){
            # my $setransexonid = $SEexonid{$tmptransid};
            # print OUTSETRANSID "$inputline\t$tmptransid\t$setransexonid\n"
            print OUTSETRANSID "$inputline\t$tmptransid\n"

        }
        
        foreach my $tmptransid(@pos1){
            
            print OUTUSTRANSID "$inputline\t$tmptransid\n"
        }

        foreach my $tmptransid(@pos2){
            
            print OUTDSTRANSID "$inputline\t$tmptransid\n"
        }

        foreach my $tmptransid(@noSEtransid){
            
            print OUTUSDSTRANSID "$inputline\t$tmptransid\n"
        }

        #  foreach my $tmptransid(@USSEDStransid){
            
        #     print OUTUSSEDSTRANSID "$inputline\t$tmptransid\n"
        # }

        ###################################################################
        ##### Start analysis contain both of US, SE and US transcripts.

        foreach my $tmptransid(@USSEDStransid){

            our @tmpexonkeys=();
            our @tmpCDSkeys=();
            
            foreach my $tmpexonkey(@sortexonkey){   #get exon keys
                if ($tmpexonkey =~ m/$tmptransid/){
                    push ( @tmpexonkeys, $tmpexonkey)
                }
            }

            foreach my $tmpCDSkey(@sortCDSkey){
                if ($tmpCDSkey =~ m/$tmptransid/){ #get CDS keys
                    push (@tmpCDSkeys, $tmpCDSkey)
                }
            }
            
            ##Get start codon

            my $startcodon_exonnumber = $trans_CDSstart{$tmptransid};
            my $stopcodon_exonnumber = $trans_CDSend{$tmptransid};

            our $setransexonid=0;

            foreach my $SEkey(@tmpexonkeys){
                # say $SEpos;
                if ($SEkey=~m/$SEpos/){
                    # say $SEkey;

                    $setransexonid=$key2exonnumber{$SEkey};
                    # say $setransexonid;

                }
            }

            # my $setransexonid = $SEexonid{$tmptransid};
            my $transexonnumbers = $trans_exonnumbers{$tmptransid};
            my $transPM = $trans_PM{$tmptransid};

            if ($transPM eq "+"){

                if ($setransexonid < $startcodon_exonnumber){

                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\t5UTR\n";

                }elsif($setransexonid > $startcodon_exonnumber and $setransexonid<$stopcodon_exonnumber){ ######### inner exons start remove SE

                    
                    our $originexonseqs = "";

                    for (our $x = 1; $x <= $transexonnumbers; $x++){
                        # say $x;
                        # say $exonseqs{$tmptransid}[$x];
                         $originexonseqs .= $exonseqs{$tmptransid}[$x];

                    }

                    our $originexonseqsLen = length($originexonseqs);

                    our $removeSEexonseqs="";

                    for (our $y = 1; $y < $setransexonid; $y++){
                        # say $y;
                        # say $exonseqs{$tmptransid}[$y];
                        $removeSEexonseqs .= $exonseqs{$tmptransid}[$y];
                    }

                    my $outSEseq = $exonseqs{$tmptransid}[$setransexonid];
                    my $outSEseqlen = length($outSEseq);

                    for (our $z = $setransexonid + 1 ; $z <= $transexonnumbers; $z++){
                        # say $z;
                        # say $exonseqs{$tmptransid}[$z];
                        $removeSEexonseqs .= $exonseqs{$tmptransid}[$z];
                    }

                    my $removeSEexonseqsLen = length($removeSEexonseqs); 
                    # Loaded exons.

                    #Start load CDS + exons.

                    our $originCDSseqs = "";

                    for (our $x = $startcodon_exonnumber; $x <= $stopcodon_exonnumber; $x++){
                        #  say $x;
                        #  say $CDSseqs{$tmptransid}[$x];
                          $originCDSseqs .= $CDSseqs{$tmptransid}[$x];

                    }

                    my $originCDSseqsLen = length($originCDSseqs);

                    our $removeSECDSseqs = "";

                    for (our $y = $startcodon_exonnumber; $y < $setransexonid; $y++){
                        # say $y;
                        # say $exonseqs{$tmptransid}[$y];
                         $removeSECDSseqs .= $CDSseqs{$tmptransid}[$y];
                         $accCDSlength{$tmptransid}[$y] = length($removeSECDSseqs);
                    }

                    my $removeSECDSseqs_upstreamLen = length($removeSECDSseqs);
                    # my $outSEseqlen = length($outSEseq);

                    for (our $z = $setransexonid + 1 ; $z <= $transexonnumbers; $z++){
                        # say $z;
                        # say $exonseqs{$tmptransid}[$z];
                         $removeSECDSseqs .= $exonseqs{$tmptransid}[$z];
                         $accCDSlength{$tmptransid}[$z] = length($removeSECDSseqs);
                    }

                    my $removeSECDSseqsLen = length($removeSECDSseqs); 

                    my $originCDSseqsAA = seq2aa1st($originCDSseqs);
                    my $removeSEexonseqsAA = seq2aa1st($removeSECDSseqs);

                    my $originCDSseqsAALen = length($originCDSseqsAA);

                    my @removeSEexonseqsAA_Pos=();


                    while ($removeSEexonseqsAA =~ m/_/gc){

                        my $tmppos = pos($removeSEexonseqsAA);
                        push (@removeSEexonseqsAA_Pos,$tmppos);

                    }

                    my $removeSEexonseqsAA_1stPos = $removeSEexonseqsAA_Pos[0];
                    my $removeSEexonseqsAA_allPos = join(', ', @removeSEexonseqsAA_Pos);

                    if (!defined($removeSEexonseqsAA_1stPos)){
                        $removeSEexonseqsAA_1stPos= "Null";
                        $removeSEexonseqsAA_allPos = "Null";

                    }                    

                    
                    $originCDSseqsAA .= "_"; #originAA add an Stop_codon.

                    my @originCDSseqsAA_Pos=();

                    while ($originCDSseqsAA =~ m/_/gc){

                        my $tmppos = pos($originCDSseqsAA);
                        push (@originCDSseqsAA_Pos,$tmppos);

                    }

                    my $originCDS_1stPos = $originCDSseqsAA_Pos[0];
                    
                    my $originCDS_allPos_allPos = join(', ', @originCDSseqsAA_Pos);

                    my $flag1 = "";

                    if($removeSEexonseqsAA_1stPos eq "Null"){
                        
                        $flag1 = "No stop_codon";

                    }elsif ($originCDS_1stPos == $removeSEexonseqsAA_1stPos){

                        $flag1 = "Need check";
                    
                    }elsif(($originCDS_1stPos * 3 - $outSEseqlen ) == ($removeSEexonseqsAA_1stPos * 3)){

                        $flag1 = "Same stop_codon";

                    }elsif(($originCDS_1stPos * 3 - $outSEseqlen ) < ($removeSEexonseqsAA_1stPos * 3)){

                        $flag1 = "Downstream stop_codon";
                        
                    }elsif(($originCDS_1stPos * 3 - $outSEseqlen ) > ($removeSEexonseqsAA_1stPos * 3)){

                        $flag1 = "Upstream stop_codon";

                    }

                    our $final_1st_stopcodon_exon_number = "-"; ######## start tell NMD.
                    our $final_1ststopcodon_accCDSlen = "-";
                    our $flag3 = "-";
                    our $NewstopcodonPos= "-";
                    our $tmpdj = "-";
                    our $final_1st_stopcodon_exon_number_len = "-";

                    if ($flag1 eq "Upstream stop_codon"){

                         $NewstopcodonPos = $removeSEexonseqsAA_1stPos * 3 ;
                        

                        for (my $a = $startcodon_exonnumber; $a <= $transexonnumbers; $a++){
                            my $tmplen; 
                            if ($a == $setransexonid){
                                next;
                            } 

                            if ($a == $startcodon_exonnumber){
                                $tmplen = $CDSsLength{$tmptransid}[$a];
                            }else{
                                $tmplen = $exonsLength{$tmptransid}[$a];
                            }

                            if (($NewstopcodonPos >= $accCDSlength{$tmptransid}[$a] - $tmplen) and ($NewstopcodonPos <= $accCDSlength{$tmptransid}[$a])){ # search 1st stopcodon length range.

                                $final_1st_stopcodon_exon_number_len = $tmplen;
                                $final_1st_stopcodon_exon_number = $a;
                                $tmpdj = $accCDSlength{$tmptransid}[$a] - $NewstopcodonPos;
                                $final_1ststopcodon_accCDSlen = $accCDSlength{$tmptransid}[$a];

                                if ($final_1st_stopcodon_exon_number == $transexonnumbers){
                                    $flag3 = "Last exon";
                                } elsif($tmpdj>= $dj){

                                    $flag3 = "NMD";  ############### NMD dj <= 50 rules.


                                }

                            }
                        }


                    }elsif ($flag1 eq "Downstream stop_codon"){

                        $NewstopcodonPos = $removeSEexonseqsAA_1stPos * 3 ;


                    }

                    


                    #End loaded CDS + exons.
                    
                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\tinner_exon";
                    print OUTUSSEDSTRANSID "\t$outSEseq\t$originexonseqs\t$removeSEexonseqs\t$outSEseqlen\t$originexonseqsLen\t$removeSEexonseqsLen\t$originCDSseqs\t$removeSECDSseqs\t$originCDSseqsLen\t$removeSECDSseqsLen";
                    print OUTUSSEDSTRANSID "\t$originCDSseqsAA\t$removeSEexonseqsAA\t$originCDSseqsAALen\t$originCDS_1stPos\t$originCDS_allPos_allPos\t$removeSEexonseqsAA_1stPos\t$removeSEexonseqsAA_allPos";
                    print OUTUSSEDSTRANSID "\t$flag1\t$final_1st_stopcodon_exon_number\t$NewstopcodonPos\t$final_1ststopcodon_accCDSlen\t$final_1st_stopcodon_exon_number_len\t$tmpdj\t$flag3\n";

                }elsif($setransexonid == $startcodon_exonnumber){

                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\tstart_codon\n";

                }elsif($setransexonid > $stopcodon_exonnumber){
                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\t3UTR\n";
                }elsif($setransexonid == $stopcodon_exonnumber){
                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\tstop_codon\n";
                }

            }elsif($transPM eq "-"){


                if ($setransexonid > $stopcodon_exonnumber){

                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\t5UTR\n";

                }elsif($setransexonid>$startcodon_exonnumber and $setransexonid<$stopcodon_exonnumber){
                
                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\tinner_exon\n";


                }elsif($setransexonid == $stopcodon_exonnumber){

                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\tstart_codon\n";

                }elsif($setransexonid < $startcodon_exonnumber){

                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\t3UTR\n";
                }elsif($setransexonid == $startcodon_exonnumber){

                    print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$transPM\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\tstop_codon\n";
                }



            }    



            # print OUTUSSEDSTRANSID "$inputline\t$tmptransid\t$setransexonid\t$transexonnumbers\t$startcodon_exonnumber\t$stopcodon_exonnumber\n";


        }


}

close INPUTLIST;
close OUTTRANSID;
close OUTSETRANSID;
close OUTUSTRANSID;
close OUTDSTRANSID;
close OUTUSDSTRANSID;
close OUTUSSEDSTRANSID;
######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;
