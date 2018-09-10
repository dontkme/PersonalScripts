#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# FindSSR v1.1001 2018/09/11
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
 our $nt=8;
 our $re=3;
#our $seqfilename ='';
#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"u=i"=>\$upstreml,"d=i"=>\$downstreml)
or die("Error in command line arguments\nUsage: perl FindSSR.pl [-o outfileprefix] <inputFASTA>
#AUTHORS
# Kaining Hu (c) 2018
# FindSSR v1.0000 2018/09/10
# hukaining\@gmail.com
\n");
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
$opfn="findSSR_out";
print "Output files: $opfn.53.txt $opfn.ssr.fa $opfn.gff3\n";
}else{
print "Output files: $opfn.53.txt $opfn.ssr.fa $opfn.gff3\n";
}


our $starttime=time();
#our $hairpinpattern="cccgccc";
#say our $testseq='((GA){4,20})';
say our $testseq='([atcg]{10}(([atcg]{1,8})(\3{3,50}))[atcg]{10})';
#say our $TCseq='([atgcn]{5}(TC(TCTACTA|T.TACTA.T|T.TACTAC|.{2}TACTACT|T.TAC.ACT|T.TA.TACT|T.TACT.CT|T.T.CTACT|.CTACTA.T|.{9}TATTAAG))[atgcn]{20})';
print "Running. Please wait for a minite.\n";
#####################################
#Start main 
#####################################
open RSUT,"> $opfn.53.txt" or die ("[-] Error: Can't open or creat $opfn.53.txt\n");
open FASTARESULT,"> $opfn.ssr.fa" or die ("[-] Error: Can't open or creat $opfn.ssr.fa\n");
#open TCRSUT, "> $opfn.5.txt" or die ("[-] Error: Can't open or creat $opfn.5.txt\n");
#open TCFARESULT, "> $opfn.5.fa" or die ("[-] Error: Can't open or creat $opfn.5.fa\n");
open OUTGFF, "> $opfn.gff3" or die ("[-] Error: Can't open or creat $opfn.gff3\n");
print OUTGFF "##gff-version 3\n";
# open DOWNFA,"> $opfn.down$downstreml.fa" or die ("[-] Error: Can't open or creat $opfn.down$downstreml.fa\n");
our $ssrcount=0;
for (our $ni=0;$ni<$Chri;$ni++){
	say "Deal with $Chrname[$ni].";
	our $count1=0;
	for (our $i=1;$i<=$nt;$i++){
		our $seqall = $Chrseq[$ni];
		if ($i == 1){
			$re=12;
			# $testseq='([atcg]{10}((([atcg])[^\3]{'.$i-1.'})(\4{'.$re.',50}))[atcg]{10})';
		} elsif($i ==2) {
			$re=6;
			#  $testseq='([atcg]{10}((([atcg])[^\3]{'.($i-1).'})(\4{'.$re.',50}))[atcg]{10})';
		} elsif ($i==3) {
			$re=5;
            # say $testseq='
            # (
            #     [ATGC]{10}
            #     (
            #         ([^A]{3}|[^T]{3}|[^G]{3}|[^C]{3}|[^N]{3})
            #         (\3{'.($re-1).',50})
            #     )[atcg]{10}
            # )';
		} else {
			$re=3;
		}
		 $testseq='([atcg]{10}(([atcg]{'.($i).'})(\3{'.($re-1).',50}))[atcg]{10})';
    ################Start deal with plus seqs.################
    while ($seqall=~ /$testseq/igcx){
 
	#say "Chrom:",$AN;
	
	#  $count1++;
	#say "No.",$count1;
	#our $l3=length ($3);
	# say "s1: $1";
    # say "s2: $2";
    # say "s3: $3";
    # say "s4: $4";
    # say "s5: $5";
    our $motif=$3;
    our $motifl=length($3);
    our $ssr=$2;
    our $ssrl=length($2);
    our $repeatn= $ssrl/$motifl;
	our $l1=length ($1);
	our $seqpos=pos($seqall);
	our $seqstarpos=$seqpos-$l1+1;
	# our $seqdownpos=$seqpos+$downstreml;
	our $ChrID=$Chrname[$ni];
    if ($ssrl>100){
        next;
	} elsif($i>=2 and $motif =~ /[A]{$i}|[T]{$i}|[G]{$i}|[C]{$i}/i){
        next;
    } elsif ($i>=4 and $motif =~ /\A([ATGC]{2,4})\1{1,3}\Z/i){
        next;
    }
    # } elsif (($motifl==1) and ($ssrl <12) ){
    #     next;
    # } elsif (($motifl==2) and ($ssrl <12) ){
    #     next;
    # } elsif (($motifl==3) and ($ssrl <15)) {
    #     next;
    # }
    $count1++;
    $ssrcount++;
	our $HID ="$ChrID"."SSR$count1";
	#say $HID;
	print RSUT  ">$HID"," ";
	print RSUT "$ChrID:$seqstarpos..$seqpos  $motif $motifl $repeatn $ssr $ssrl\n";
	print RSUT "$1\n";
 #	print RSUT "-----------------------------------------\n";
	#print FASTARESULT ">$ChrID"."H$count1.up $ChrID:$seqstarpos..$seqpos\n";
	# our $rightfinalseq="";
	# if ($seqpos<$upstreml){
	#  $rightfinalseq = substr($seqall,1,$seqpos);
	# }else{
	#  $rightfinalseq = substr($seqall,$seqpos-$upstreml,$upstreml);
	# }
	# print FASTARESULT "$rightfinalseq\n";
	
	# print DOWNFA ">$ChrID"."H$count1.down $ChrID:$seqpos..$seqdownpos\n";
	# our $downfasta = substr($seqall,$seqpos,$downstreml);
	# print DOWNFA "$downfasta\n";
	
	print OUTGFF "$ChrID\tFindSSR\tCDS\t".($seqstarpos+10)."\t".($seqpos-10)."\t.\t+\t.\tName=$HID;ID=$HID;Parent=$HID;Type=$motif;Super_Family=$repeatn\n";
	print FASTARESULT ">$HID $ChrID:".($seqstarpos+10)."..".($seqpos-10)."  $motif $motifl $repeatn $ssr $ssrl\n";
	print FASTARESULT "$ssr\n";
	# our $TCi=0;
	# #say our $TCseq='([atgcn]{10}(TC(TCTACTA|T.TACTA.T|T.TACTAC|.{2}TACTACT|T.TAC.ACT|T.TA.TACT|T.TACT.CT|T.T.CTACT|.CTACTA.T|.{9}TATTAAG))[atgcn]{15})';
	# our $rfsl=length($rightfinalseq);
	
	
    }
	}
    ################ deal with plus seqs End#################




# ################Start TR-seqs#############################
# 	 $seqall = TRseq($seqall);
# 	 our $count2=0;
# 	our $seqalll = length($seqall);
# while ($seqall=~ /$testseq/igc){
# #while ($seqall =~ /$testseq/igc){
# 	#say $seqfilename;
# 	#say "Chrom:",$AN;
# 	 $count2++;
# 	#say "No.",$count2;
# 	#our $l3=length ($3);
# 	our $l1=length ($1);
# 	our $seqpos=$seqalll-(pos($seqall))+1;
# 	our $seqstarpos=$seqpos+$l1-1;
# 	our $ChrID=$Chrname[$ni];

# 	our $trHID ="tr$ChrID"."H$count2";
# 	#say $trHID;
# 	print RSUT  ">tr$ChrID","H",$count2," ";
# 	print RSUT "$ChrID:$seqstarpos..",$seqpos," \n";
# 	#print RSUT $l3," ";
# 	#print RSUT $2," ";
# 	#print  RSUT $3," ";
# 	#print  RSUT $4,"\n";
# 		print RSUT $1,"\n";
# #	print RSUT "-----------------------------------------\n";
# 	print FASTARESULT ">tr$ChrID"."H$count2.up $ChrID:$seqstarpos..$seqpos\n";
# 	our $trseqpos = pos($seqall);
# 	our $rightfinalseq="";
# 	if($trseqpos<$upstreml){
# 	 $rightfinalseq = substr($seqall,1,$trseqpos);
# 	}else{
# 	 $rightfinalseq = substr($seqall,$trseqpos-$upstreml,$upstreml);
# 	}
# 	print FASTARESULT "$rightfinalseq\n";
# 	if ($seqpos<$downstreml) {
#     print DOWNFA ">tr$ChrID"."H$count2.down $ChrID:$seqpos..1\n";
#   }else {
#     my $temppos=$seqpos-$downstreml;
#     print DOWNFA ">tr$ChrID"."H$count2.down $ChrID:$seqpos..$temppos\n";
#   }
# 	#print DOWNFA ">tr$ChrID"."H$count2.down $ChrID:$seqpos..$seqpos-$downstreml\n";
# 	our $downfasta = substr($seqall,$trseqpos,$downstreml);
# 	print DOWNFA "$downfasta\n";
	
	
# 	print OUTGFF "$ChrID\tFindPattern\tCDS\t$seqpos\t$seqstarpos\t.\t-\t.\tName=$trHID.3;ID=$trHID;Parent=$trHID;Type=GA;Super_Family=GA\n";

# 	our $TCi=0;
# 	#say our $TCseq='([atgcn]{10}(TC(TCTACTA|T.TACTA.T|T.TACTAC|.{2}TACTACT|T.TAC.ACT|T.TA.TACT|T.TACT.CT|T.T.CTACT|.CTACTA.T|.{9}TATTAAG))[atgcn]{15})';
# 	our $rfsl=length($rightfinalseq);

	
# }	
# ################ TR-seqs End#############################

}   


close RSUT;
#close TCRSUT; 
close FASTARESULT;
#close TCFARESULT;
close OUTGFF;
# close DOWNFA;

######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! Find $ssrcount SSR record(s). Used %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;