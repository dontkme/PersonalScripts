#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Find Pattern (EAHelitron)v1.0000 2018/07/26
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
our $opfn="";
my $verbose;
our $upstreml=3000;
our $downstreml=1100;
#our $seqfilename ='';
#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"u=i"=>\$upstreml,"d=i"=>\$downstreml)
or die("Error in command line arguments\nUsage: perl FindPattern [-o outfileprefix][-u|upstremlength int][-d|downstremlength int] <inputFASTA>\n");
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

our $loadingstarttime=time();

print "Start loading genomeic sequence.\n";
#if(!open SEQFILENAME,"< $seqfilename"){
# die "File not Found\n";
#}
#if(!open SEQFILENAME,"< $_"){
# die "File not Found\n";
#}
our $Chri=0;
our @Chrname=();
our @Chrseq=();
#@ARGV = qw#''  Not_Find_a_File#;
#say @ARGV;
#say $0;
while(defined(our $seq = <>)){
#while(our $seq = <SEQFILENAME>){
	if ($seq =~ m/^.*>/) {
	$seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
	print "$1\n";
	 $Chrname[$Chri]= $1;
	$Chri++;
	}else{
		$seq =~ s/\s//;
		$seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
	 $Chrseq[$Chri-1] .=$seq;
			}
}
#close SEQFILENAME;
our $loadingendtime=time();
print "$Chri Sequences\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;
if ($opfn eq ""){
$opfn="EAHeli_out";
print "Output files:$opfn.fas $opfn.gff3\n";
}else{
print "Output files:$opfn.fas $opfn.gff3\n";
}


our $starttime=time();
our $hairpinpattern="cccgccc";
say our $testseq='((GA){4,20})';
#say our $TCseq='([atgcn]{5}(TC(TCTACTA|T.TACTA.T|T.TACTAC|.{2}TACTACT|T.TAC.ACT|T.TA.TACT|T.TACT.CT|T.T.CTACT|.CTACTA.T|.{9}TATTAAG))[atgcn]{20})';
print "Running. Please wait for a minite.\n";
#####################################
#Start main 
#####################################
open RSUT,"> $opfn.3.txt" or die ("[-] Error: Can't open or creat $opfn.3.txt\n");
open FASTARESULT,"> $opfn.u$upstreml.fa" or die ("[-] Error: Can't open or creat $opfn.u$upstreml.fa\n");
#open TCRSUT, "> $opfn.5.txt" or die ("[-] Error: Can't open or creat $opfn.5.txt\n");
#open TCFARESULT, "> $opfn.5.fa" or die ("[-] Error: Can't open or creat $opfn.5.fa\n");
open OUTGFF, "> $opfn.gff3" or die ("[-] Error: Can't open or creat $opfn.gff3\n");
print OUTGFF "##gff-version 3\n";
open DOWNFA,"> $opfn.down$downstreml.fa" or die ("[-] Error: Can't open or creat $opfn.down$downstreml.fa\n");

for (our $ni=0;$ni<$Chri;$ni++){
	say "Deal with $Chrname[$ni].";
	our $count1=0;
	our $seqall = $Chrseq[$ni];
################Start deal with plus seqs.################
while ($seqall=~ /$testseq/igc){
#while ($seqall =~ /$testseq/igc){
	#say $seqfilename;
	#say "Chrom:",$AN;
	 $count1++;
	#say "No.",$count1;
	#our $l3=length ($3);
	our $l1=length ($1);
	our $seqpos=pos($seqall);
	our $seqstarpos=$seqpos-$l1+1;
	our $seqdownpos=$seqpos+$downstreml;
	our $ChrID=$Chrname[$ni];

	our $HID ="$ChrID"."H$count1";
	#say $HID;
	print RSUT  ">$HID"," ";
	print RSUT "$ChrID:$seqstarpos..",$seqpos," \n";
	#print RSUT $l3," ";
	#print RSUT $2," ";
	#print  RSUT $3," ";
	#print  RSUT $4,"\n";
	print RSUT $1,"\n";
#	print RSUT "-----------------------------------------\n";
	print FASTARESULT ">$ChrID"."H$count1.up $ChrID:$seqstarpos..$seqpos\n";
	our $rightfinalseq="";
	if ($seqpos<$upstreml){
	 $rightfinalseq = substr($seqall,1,$seqpos);
	}else{
	 $rightfinalseq = substr($seqall,$seqpos-$upstreml,$upstreml);
	}
	print FASTARESULT "$rightfinalseq\n";
	
	print DOWNFA ">$ChrID"."H$count1.down $ChrID:$seqpos..$seqdownpos\n";
	our $downfasta = substr($seqall,$seqpos,$downstreml);
	print DOWNFA "$downfasta\n";
	
	print OUTGFF "$ChrID\tFindpattern\tCDS\t$seqstarpos\t$seqpos\t.\t+\t.\tName=$HID.3;ID=$HID;Parent=$HID;Type=GA;Super_Family=GA\n";
	
	our $TCi=0;
	#say our $TCseq='([atgcn]{10}(TC(TCTACTA|T.TACTA.T|T.TACTAC|.{2}TACTACT|T.TAC.ACT|T.TA.TACT|T.TACT.CT|T.T.CTACT|.CTACTA.T|.{9}TATTAAG))[atgcn]{15})';
	our $rfsl=length($rightfinalseq);
	
	
}
################ deal with plus seqs End#################




################Start TR-seqs#############################
	 $seqall = TRseq($seqall);
	 our $count2=0;
	our $seqalll = length($seqall);
while ($seqall=~ /$testseq/igc){
#while ($seqall =~ /$testseq/igc){
	#say $seqfilename;
	#say "Chrom:",$AN;
	 $count2++;
	#say "No.",$count2;
	#our $l3=length ($3);
	our $l1=length ($1);
	our $seqpos=$seqalll-(pos($seqall))+1;
	our $seqstarpos=$seqpos+$l1-1;
	our $ChrID=$Chrname[$ni];

	our $trHID ="tr$ChrID"."H$count2";
	#say $trHID;
	print RSUT  ">tr$ChrID","H",$count2," ";
	print RSUT "$ChrID:$seqstarpos..",$seqpos," \n";
	#print RSUT $l3," ";
	#print RSUT $2," ";
	#print  RSUT $3," ";
	#print  RSUT $4,"\n";
		print RSUT $1,"\n";
#	print RSUT "-----------------------------------------\n";
	print FASTARESULT ">tr$ChrID"."H$count2.up $ChrID:$seqstarpos..$seqpos\n";
	our $trseqpos = pos($seqall);
	our $rightfinalseq="";
	if($trseqpos<$upstreml){
	 $rightfinalseq = substr($seqall,1,$trseqpos);
	}else{
	 $rightfinalseq = substr($seqall,$trseqpos-$upstreml,$upstreml);
	}
	print FASTARESULT "$rightfinalseq\n";
	if ($seqpos<$downstreml) {
    print DOWNFA ">tr$ChrID"."H$count2.down $ChrID:$seqpos..1\n";
  }else {
    my $temppos=$seqpos-$downstreml;
    print DOWNFA ">tr$ChrID"."H$count2.down $ChrID:$seqpos..$temppos\n";
  }
	#print DOWNFA ">tr$ChrID"."H$count2.down $ChrID:$seqpos..$seqpos-$downstreml\n";
	our $downfasta = substr($seqall,$trseqpos,$downstreml);
	print DOWNFA "$downfasta\n";
	
	
	print OUTGFF "$ChrID\tFindPattern\tCDS\t$seqpos\t$seqstarpos\t.\t-\t.\tName=$trHID.3;ID=$trHID;Parent=$trHID;Type=GA;Super_Family=GA\n";

	our $TCi=0;
	#say our $TCseq='([atgcn]{10}(TC(TCTACTA|T.TACTA.T|T.TACTAC|.{2}TACTACT|T.TAC.ACT|T.TA.TACT|T.TACT.CT|T.T.CTACT|.CTACTA.T|.{9}TATTAAG))[atgcn]{15})';
	our $rfsl=length($rightfinalseq);

	
}	
################ TR-seqs End#############################
}
close RSUT;
#close TCRSUT; 
close FASTARESULT;
#close TCFARESULT;
close OUTGFF;
close DOWNFA;

######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;