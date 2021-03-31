#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2021
# Get Sequences from GTF with genome v1.0000 2021/03/29
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
our $opfn="getseqsOut";
my $verbose;
our $upstreml=5000;
our $downstreml=5000;
our $annot;
our $sortID="gene_id";
our $feature="CDS";
our $inputlist="./input.txt";

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"u=i"=>\$upstreml,"d=i"=>\$downstreml, "g=s"=>\$annot,"f=s"=>\$feature,"s=s"=>\$sortID, "i=s"=>\$inputlist)
or die("[-]Error in command line arguments
  Usage: perl GetseqsfromGTF [options] <-g string|GTF annoation file> <input FASTA file(s)>
    options:
    [-o string|outprefix Default: getseqsOut]
    [-s string|Specify attribute type in GFF annotation for sorting. default: gene_id]
    [-f string|Specify feature type in GFF annotation.default: CDS]
    [-u int|upstream length Default: 5000]
    [-d int|downstream length Default: 5000]
    [-i input list]
	 
    Note: Get Sequences from GTF with genome v1.0000 2021/03/29.\n");

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

open OUT, "> $opfn.fa" or die ("[-] Error: Can't open or creat $opfn.fa\n");

open INPUTLIST, "< $inputlist" or die ("[-] Error: Can't open list file: $inputlist.\n");

open OUTSEQ, "> $opfn.outseq.txt" or die ("[-] Error: Can't open or creat $opfn.outseq.txt\n");

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
print "Start loading input GTF.\n";
print "ID: $sortID\n";
print "Feature: $feature\n";


    

    our @tmp;
    our $faheadid="test";
    our $annotcount=0;
    our (%key2Chr,%key2type, %key2startpos, %key2endpos, %key2PM,%key2restwords,%key2line);

while(defined(our $inrow = <ANNOT>)){

    if ($inrow =~ m/^\#/) {next;}
    if ($annotcount > 1 and $annotcount % 1000 == 0){
        print "Dealed with $annotcount annotations.\n";
    }
    @tmp = split (/\t/,$inrow);
    # say $inrow;
    
    # if ($tmp[2] ne $feature){next;}
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
        }elsif($tmp3=~m/transcript_id \"(\S+)\"/i){
            $transcript_id = $1;            
        }elsif($tmp3 =~ m/exon_number \"(\S+)\"/i){
            $exon_number = $1 + 0;
            # $exon_numberf3 = printf "%3.0f",$exon_number;
            $exon_numberf3 = sprintf "%03d",$exon_number;
        }elsif($tmp3 =~ m/exon_id \"(\S+)\"/i){
            $exon_id = $1;
        }elsif($tmp3 =~ m/gene_name \"(\S+)\"/i){
            $gene_name = $1;
        }
    }

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



    # if ($plusminus eq "+"){
    #     my $finalseq="";
    #     $finalseq= substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1);
    #     print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
    #     print OUT "$finalseq\n";

    #     my $upfinalseq="";
    #     if ($seqstartpos<$upstreml){
    #         print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:1..$seqstartpos $plusminus\n";
    #         $upfinalseq = substr($Chrid2seq{$seqchrid},0,$seqstartpos);
    #     }else{
    #         print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:".($seqstartpos-$upstreml)."..$seqstartpos $plusminus\n";
    #         $upfinalseq = substr($Chrid2seq{$seqchrid},$seqstartpos-$upstreml,$upstreml);

    #     }
       
    #     print OUTUP "$upfinalseq\n";

    #     my $downfinalseq="";
    #     $downfinalseq = substr($Chrid2seq{$seqchrid},$seqendpos,$downstreml);
    #     print OUTDOWN ">$faheadid.down$downstreml"."  $seqchrid:$seqendpos..".($seqstartpos+$downstreml)." $plusminus\n";
    #     print OUTDOWN "$downfinalseq\n";
        
    #     my $useqd="";
    #     if ($seqstartpos<$upstreml){
    #         print OUTALL ">$faheadid.u$upstreml.d$downstreml"."  $seqchrid:1..".($seqendpos+$downstreml)." $plusminus\n";
    #         $useqd = substr($Chrid2seq{$seqchrid},0,($seqendpos+$downstreml));
    #     }else{
    #         print OUTALL ">$faheadid.u$upstreml.d$downstreml"."  $seqchrid:".($seqstartpos-$upstreml)."..".($seqendpos+$downstreml)." $plusminus\n";
    #         $useqd = substr($Chrid2seq{$seqchrid},($seqstartpos-$upstreml-1),($upstreml+$seqendpos-$seqstartpos+1+$downstreml));

    #     }
    #     print OUTALL "$useqd\n";
    
    # } elsif($plusminus eq "-"){
    #     my $finalseq= TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1));
    #     print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
    #     print OUT "$finalseq\n";

    #     my $downfinalseq="";
    #     if ($seqstartpos<$downstreml){
    #         print OUTDOWN ">$faheadid.down$downstreml"."  $seqchrid:1..$seqstartpos $plusminus\n";
    #         $downfinalseq = TRseq(substr($Chrid2seq{$seqchrid},0,$seqstartpos-1));
    #         print OUTDOWN "$downfinalseq\n";
    #     }else{
    #         print OUTDOWN ">$faheadid.down$downstreml"."  $seqchrid:".($seqstartpos-$downstreml)."..$seqstartpos $plusminus\n";
    #         $downfinalseq = TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos-$downstreml-1,$downstreml));
    #         print OUTDOWN "$downfinalseq\n";

    #     }
    #     # print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:$seqstartpos..$seqendpos\n";
        

    #     my $upfinalseq="";
    #     $upfinalseq = TRseq(substr($Chrid2seq{$seqchrid},$seqendpos,$upstreml));
    #     print OUTUP ">$faheadid.up$upstreml"."  $seqchrid:$seqendpos..".($seqendpos+$upstreml)." $plusminus\n";
    #     print OUTUP "$upfinalseq\n";

    #     my $useqd="";
    #     if ($seqstartpos<$downstreml){
    #         print OUTALL ">$faheadid.u$upstreml.d$downstreml"."  $seqchrid:1..".($seqendpos+$upstreml)." $plusminus\n";
    #         $useqd = TRseq(substr($Chrid2seq{$seqchrid},0,($seqendpos+$upstreml-1)));
    #         print OUTALL "$useqd\n";
    #     }else{
    #         print OUTALL ">$faheadid.u$upstreml.d$downstreml"."  $seqchrid:".($seqstartpos-$downstreml)."..".($seqendpos+$upstreml)." $plusminus\n";
    #         $useqd = TRseq(substr($Chrid2seq{$seqchrid},($seqstartpos-$downstreml-1),($upstreml+$seqendpos-$seqstartpos+1+$downstreml)));
    #         print OUTALL "$useqd\n";

    #     }

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

our @sortkey = (sort keys %key2PM);

foreach my $sortkey(@sortkey){
    # say $sortkey;
    # if ($sortkey =~m/exon/i){
    if ($sortkey =~m/$feature/i){
        my $finalseq="";
        my $seqchrid=$key2Chr{$sortkey};
        my $seqstartpos=$key2startpos{$sortkey};
        my $seqendpos=$key2endpos{$sortkey};
        my $PM=$key2PM{$sortkey};
        $finalseq = getseq($sortkey,$seqchrid,$seqstartpos,$seqendpos,$PM);
        print OUT ">$sortkey\n";
        print OUT "$finalseq\n";
    #  if ($sortkey =~m//i){
    #     if ($key2PM{$sortkey} eq "+"){
    #         my $finalseq="";
    #         my $seqchrid=$key2Chr{$sortkey};
    #         my $seqstartpos=$key2startpos{$sortkey};
    #         my $seqendpos=$key2endpos{$sortkey};
            
    #         $finalseq= substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1);
    #         # print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
    #         print OUT ">$sortkey\n";
    #         print OUT "$finalseq\n";

    #     }elsif($key2PM{$sortkey} eq "-"){
    #             my $seqchrid=$key2Chr{$sortkey};
    #             my $seqstartpos=$key2startpos{$sortkey};
    #             my $seqendpos=$key2endpos{$sortkey};
    #             my $finalseq="";
    #             $finalseq= TRseq(substr($Chrid2seq{$seqchrid},$seqstartpos-1,$seqendpos-$seqstartpos+1));
    #             # print OUT ">$faheadid"." $seqchrid:$seqstartpos..$seqendpos $plusminus\n";
    #             print OUT ">$sortkey\n";
    #             print OUT "$finalseq\n";

    #     }
    }
}
# say @sortkey;







close OUT;

######################## condon2aa #######################

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

####################### codon2aa 



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
# close OUTDOWN;
# close OUTUP;

######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;
