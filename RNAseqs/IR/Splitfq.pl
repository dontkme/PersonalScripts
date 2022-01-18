#!/bin/perl
use strict;
use 5.010;
#use Getopt::Long;
#use File::Basename;
#use lib dirname $0;
#open IN, $ARGV[0];
our $count=0;
while(<>){
        $count++;
        chomp(my $readName=$_);
        #my @temp=split “:”, $readName;
        chomp (my $seqStr=<>);
        chomp (my $name2=<>);
        chomp (my $qualStr=<>);

        if($count%2 == 1){
                say $readName;
                say $seqStr;
                say $name2;
                say $qualStr;
        }else{
                say STDERR $readName;
                say STDERR $seqStr;
                say STDERR $name2;
                say STDERR $qualStr;
        }
}
