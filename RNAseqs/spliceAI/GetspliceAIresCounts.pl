#!/usr/bin/perl

use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';

our $opfn="";
our $vcfgz="";
our $suffix="res";
GetOptions("o=s" => \$opfn, "g=s" => \$vcfgz,"s=s" => \$suffix)
or die("[-]Error in command line arguments
  Usage: perl GetspliceAIresCounts.pl [options] -g <in.vcf.gz> <input region file>
    options:
	 [-o string|output prefix Default: Res_Counts]
     [-s string|output suffix Default: res]
    Note:  Get spliceai results count v0.003 2022/02/07.\n");


if ($opfn eq ""){
	$opfn="Res_Counts";
	print "output prefix:$opfn\n";
}else{
	print "output prefix:$opfn\n";
}

if ($vcfgz eq ""){
    die ("[-]Error. Not find a vcf.gz file.\n");
}

### Main
our $starttime=time();
open OUTSUM, "> $opfn.txt" or die ("[-] Error: Can't open or create $opfn.txt\n");

while(defined(our $seq = <>)){
  if ($seq =~ m/^#.*/) {
	
    next;
  }else {
    my @chrr=split(/\t/,$seq);
    my $gene=$chrr[0];
    my $chr=$chrr[1];
    $chr=~s/chr//; ### for diff chr header
    my $pos1=$chrr[2];
    my $pos2=$chrr[3];
    $pos2 =~ s/\s//g;
    my $AS=join("_",$gene,$chr,$pos1,$pos2);
    my $ASregion="$chr:$pos1-$pos2";
    
    system "mkdir $AS";
    # system "cd $AS";
    # system "bcftools query -H -f \'\%CHROM\t\%POS\t\%ID\t\%REF\t\%ALT\t\%QUAL\t\%FILTER\t\%INFO/AC\t\%INFO/AN\t\%INFO/AF\t\%INFO/SpliceAI\n\' -r $ASregion $vcfgz > $AS/$AS.res.txt";
    # system "bcftools query -H -f \'\%CHROM\t\%POS\t\%ID\t\%REF\t\%ALT\t\%QUAL\t\%FILTER\t\%INFO/AC\t\%INFO/AN\t\%INFO/AF\n\' -r $ASregion $vcfgz > $AS/$AS.$suffix.txt";
    system "bcftools query -H -f \'\%CHROM\t\%POS\t\%ID\t\%REF\t\%ALT\t\%QUAL\t\%FILTER\t\%INFO\n\' -r $ASregion $vcfgz > $AS/$AS.$suffix.txt";
    # system "bcftools query -H -r $ASregion $vcfgz > $AS.res.txt";

    # my $rescount = system("wc -l $AS/$AS.res.txt | cut -c 2");
    my $rescount=0;
    open(FILE, "< $AS/$AS.$suffix.txt") or die "can't open $AS/$AS.$suffix.txt: $!";
    $rescount++ while <FILE>;
    # say $rescount;
    $rescount=$rescount-1;
    # say $rescount;
	# system 'cd ..';
    print OUTSUM "$AS\t$ASregion\t$gene\t$chr\t$pos1\t$pos2\t$rescount\n";
  }

}
close OUTSUM;
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;