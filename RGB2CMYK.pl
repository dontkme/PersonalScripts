#!/bin/perl

#AUTHORS
# Kaining Hu (c) 2017
# Convert RGB to CMYK v1.0000 2017/12/14
# hukaining@gmail.com
#
use strict;
use warnings;
#use 5.010;
use List::Util qw[min max];
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';

our $opfn="";

my $verbose;



GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
or die("Error in command line arguments\nUsage: 	perl Cutfa [options] <inputfile>\n
	 options:\n
	 [-o string|output prefix default: testup]\n
	 
	 Note: Convert RGB to CMYK v1.0000 2017/06/07\n");



#print $ARGV[0],"\n";
if (not $ARGV[0]) {
	die ("[-] Error: Not find a inputfile.\n");
}
my @tmpfliename=@ARGV;
if ($opfn eq ""){
	$opfn="$ARGV[0]2CMYK";
	print "Output file: $opfn.txt \n";
}else{
	print "Output file: $opfn.txt \n";
}

open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or creat $opfn.fa\n");

our $starttime=time();

our $tmpseqname="";
our $count1=0;
our $count2=0;


######main
print OUT "C\tM\tY\tK\n";
while(defined(our $RGBfile=<>)){
  $count1++;
  my @RGB=split/\t/,$RGBfile;
  print $RGB[0]."\n";
  my $R=$RGB[0]/255;
  my $G=$RGB[1]/255;
  my $B=$RGB[2]/255;
     my $c = 1 - $R;
    my $m = 1 - $G;
    my $y = 1 - $B;
    my $k = min($c, $m, $y);
    
    if ($k == 1) {
        $c = $m = $y = 0;   # Pure black
    } else {
        $c = ($c - $k) / (1 - $k);
        $m = ($m - $k) / (1 - $k);
        $y = ($y - $k) / (1 - $k);
    }
    print OUT "$c\t$m\t$y\t$k\n";
    $count2++;
	} #while End.


#print "Finished Cut @tmpfliename upstream $upstreml bp to $opfn.fa\nStrict mode: $strict \n$count1 seqs cut $count2 seqs.\n";
close OUT;
our $endtime=time();
printf "Done! load $count1 RGB, out $count2 CMYK.\n %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;