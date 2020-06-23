#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2020
# Combine 4 cds fasta by input files v1.0000 2020/06/23
# hukaining@gmail.com



#use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
my $verbose;

our $cds1 = "./AtCDS.txt";
our $cds2 = "./BrCDS.txt";
our $cds3 = "./BoCDS.txt";
our $cds4 = "./BnCDS.txt";
our $cds1prefix = "";
our $cds2prefix = "";
our $cds3prefix = "";
our $cds4prefix = "";

our $num=1;
our $cdsname="";
our %seqfa1;
our %seqfa2;
our %seqfa3;
our %seqfa4;
our $seqfa="seqfa";
our $opdir = "./Get4seqcdsOut";


GetOptions("o=s" => \$opdir,"verbose"=>\$verbose,"c1=s"=>\$cds1,"c2=s"=>\$cds2,"c3=s"=>\$cds3,"c4=s"=>\$cds4,"c1p=s"=>\$cds1prefix,"c2p=s"=>\$cds2prefix,"c3p=s"=>\$cds3prefix,"c4p=s"=>\$cds4prefix,)
or die("Error in command line arguments\nUsage: 	perl Get4seqcds.pl [options] <Input 4seq name list file>\n
	 options:\n
	 [-o string|output dir default: ./Get4seqcdsOut]\n
	 [-c1 string|CDS fasta file1]\n
	 [-c2 string|CDS fasta file2]\n
	 [-c3 string|CDS fasta file3]\n
	 [-c4 string|CDS fasta file4]\n
	 [-c1p string|CDS 1 prefix]\n
	 [-c2p string|CDS 2 prefix]\n
	 [-c1p string|CDS 3 prefix]\n
	 [-c1p string|CDS 4 prefix]\n

	 Note: Combine 4 cds fasta by input files v1.0000 2020/06/23\n");


our $loadingstarttime=time();

mkdir $opdir unless -d $opdir;

##########################
# lodading function
############################
sub loadingcds ($){
    my $num=$_[0];
    our $CDS=join('',"CDS",$num);
    say "$CDS";
    our $seqfa=join('',"seqfa",$num);
    say "$seqfa";
    %{$seqfa};
    our $cds=join('',"cds",$num);
    say $cds;
    our $cdsprefix=join('',"cds",$num,"prefix");
    say $cdsprefix;

    if(!open $CDS,"< ${$cds}"){
    die "${$cds} File not Found\n";
    }

    my $seq1="";
    while($seq1 = <$CDS>){
    if ($seq1 =~ m/^.*>/) {
        $seq1=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
        $cdsname = $1;
        ${$seqfa}{$cdsname}=">$$cdsprefix"."$cdsname\n";
        }else {
        ${$seqfa}{$cdsname}.=$seq1;
        }
    }

    close $CDS;
}

#############
# lodading function end
###############

for (our $i=1;$i<5;$i++){
    &loadingcds($i);
}




our $count1=0;

while(defined(our $seq = <>)){
#while(our $seq = <SEQFILENAME>){
	# my @sep = split/\t/,$seq;
	# $sep[1] =~ s/\s//g;
	# $sep[0] =~ s/\s//g;
	$seq =~ s/\s//g;
    # say $seq[0];
    say $seq;
    $sep="$seq"."_";
	# my $seqname="$sep[0]_2_"."$sep[1]";
	my $seqname=join("$sep",$cds1prefix,$cds2prefix,$cds3prefix,$cds4prefix,"");
	open OUTFA, "> $opdir/$seqname.fa" or die ("[-] Error: Can't open or creat $opdir/$seqname.fa\n");
	
  print OUTFA "$seqfa1{$seq}";
  print OUTFA "$seqfa2{$seq}";
  print OUTFA "$seqfa3{$seq}";
  print OUTFA "$seqfa4{$seq}";
  #print OUTFA "$seqfa{$sep[1]}";
  
	close OUTFA;
  $count1 ++;
			}
print "Finished combining 4 CDS! $count1 fas in $opdir\n";
our $loadingendtime=time();
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;

