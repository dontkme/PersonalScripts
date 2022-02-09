#!/usr/bin/perl

use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
use List::Util qw( min max );

our $opfn="";
our $vcfgz="";
our $suffix="res";
our $spliceai_res="NO";
our $cutoff=0.5;
GetOptions("o=s" => \$opfn, "g=s" => \$vcfgz,"s=s" => \$suffix,"r=s" => \$spliceai_res,"c=f" => \$cutoff)
or die("[-]Error in command line arguments
  Usage: perl GetspliceAIresCounts.pl [options] -g <in.vcf.gz> <input region file>
    options:
	 [-o string|output prefix Default: Res_Counts]
   [-s string|output suffix Default: res]
   [-r string|whether input vcf.gz contain spliceAI annotations. (Y)es or (N)o. Default: NO]
   [-c float|Cutoff of delta Score. Default: 0.5]
    Note:  Get spliceai results count (ClinVar version) v0.120 2022/02/09.\n");


if ($opfn eq ""){
	$opfn="Res_Counts";
	print "Output summary: $opfn.txt\n";
}else{
	print "Output summary: $opfn.txt\n";
}

if ($vcfgz eq ""){
    die ("[-]Error. Not find a vcf.gz file.\n");
}else{
  print "Input vcf.gz: $vcfgz\n";
}


if ($spliceai_res =~ m/[Y|YES]/i){
  print "Input vcf.gz contain SpliceAI annotation.\n";
}else{
  print "Input vcf.gz NOT contain SpliceAI annotation.\n";
}
say "Suffix: $suffix";
# say "Cutoff: $cutoff";
if ($cutoff > 1 && $cutoff < 0){
    die ("[-]Error. Cutoff range must from 0 to 1.\n");
}else{
  print "Cutoff: $cutoff\n";
}
### Main
our $starttime=time();
open OUTSUM, "> $opfn.txt" or die ("[-] Error: Can't open or create $opfn.txt\n");
our $inputcount=0;

while(defined(our $seq = <>)){
  if ($seq =~ m/^#.*/) {
	
    next;
  }else {
    $inputcount++;
    my @chrr=split(/\t/,$seq);
    my $gene=$chrr[0];
    my $chr=$chrr[1];
    $chr=~s/chr//; ### for diff chr header
    my $pos1=$chrr[2];
    my $pos2=$chrr[3];
    my $PM=$chrr[4];
    # $pos2 =~ s/\s//g;
    $PM =~ s/\s//g;
    my $AS=join("_",$gene,$chr,$pos1,$pos2);
    my $ASregion="$chr:$pos1-$pos2";
    
    system "mkdir -p $AS";
    # system "cd $AS";
    # system "bcftools query -H -f \'\%CHROM\t\%POS\t\%ID\t\%REF\t\%ALT\t\%QUAL\t\%FILTER\t\%INFO/AC\t\%INFO/AN\t\%INFO/AF\t\%INFO/SpliceAI\n\' -r $ASregion $vcfgz > $AS/$AS.res.txt";
    # system "bcftools query -H -f \'\%CHROM\t\%POS\t\%ID\t\%REF\t\%ALT\t\%QUAL\t\%FILTER\t\%INFO/AC\t\%INFO/AN\t\%INFO/AF\n\' -r $ASregion $vcfgz > $AS/$AS.$suffix.txt";
    if ($spliceai_res =~ m/[Y|YES]/i){

      system "bcftools query -H -f \'\%CHROM\t\%POS\t\%ID\t\%REF\t\%ALT\t\%QUAL\t\%FILTER\t\%INFO/SpliceAI\n\' -r $ASregion $vcfgz > $AS/$AS.$suffix.txt";

      
    }else {
 
      system "bcftools query -H -f \'\%CHROM\t\%POS\t\%ID\t\%REF\t\%ALT\t\%QUAL\t\%FILTER\t\%INFO\n\' -r $ASregion $vcfgz > $AS/$AS.$suffix.txt";
    
    }
    # system "bcftools query -H -r $ASregion $vcfgz > $AS.res.txt";

    # my $rescount = system("wc -l $AS/$AS.res.txt | cut -c 2");

    #### Count line numbers.
    my $rescount=0;
    my $passcount=0;
    my $maxDS=0;
    my @DS=();
    open(FILE, "< $AS/$AS.$suffix.txt") or die "can't open $AS/$AS.$suffix.txt: $!";
    open OUTPASS, "> $AS/$AS.$suffix.PASS.txt" or die ("[-] Error: Can't open or create $AS/$AS.$suffix.PASS.txt\n");


    while(defined(our $resline = <FILE>)){
        $rescount++ ;
        if ($resline =~ m/^#.*/) {
	
          next;

        }else{
          if ($spliceai_res =~ m/[Y|YES]/i){
          # if ($resline=~ m/SpliceAI=(.*)$/){
            my @sepline=split(/\t/,$resline);
            my $snppos=$sepline[1]; #### Add PASS distence calculate.
            my $spanno=$sepline[7];
            $spanno =~ s/\s//g;
            # say $spanno;
            if ($spanno eq "."){  ## Skip "." spliceAI results.
              # say $sepline[0];
              # say $sepline[1];
              next;
            }
            # say $spanno;
            my @spres=split(/\|/,$spanno);

            if ($spres[2] eq "."){  ## Skip "." spliceAI results.
             
              next;

            }
            
            my $AG=$spres[2]+0;
            my $AL=$spres[3]+0;
            my $DG=$spres[4]+0;
            my $DL=$spres[5]+0;

            my $AGP=$spres[6];
            my $ALP=$spres[7];
            my $DGP=$spres[8];
            my $DLP=$spres[9];
            $DLP =~ s/\s//g;
            # $DL =~ s/\s//g;
            @DS=($AG,$AL,$DG,$DL);
            if (scalar @DS == 0){
              next;
              }else{
              # print @DS;

              $maxDS = max(@DS);
              # say $AG;
              # say $maxDS;

              if($maxDS gt $cutoff){
                $passcount++;
                $resline =~ s/\s$//;
                my $Dpos1=$snppos-$pos1;
                my $Dpos2=$snppos-$pos2;
                print OUTPASS "$resline\t$AG\t$AL\t$DG\t$DL\t$AGP\t$ALP\t$DGP\t$DLP\t$maxDS\t$gene\t$chr\t$pos1\t$pos2\t$PM\t$Dpos1\t$Dpos2\n";
              }

            }

        # }
        }
      }
    }
    # say $rescount;
    $rescount=$rescount-1;
    # say $rescount;
	# system 'cd ..';
    print OUTSUM "$AS\t$ASregion\t$gene\t$chr\t$pos1\t$pos2\t$rescount\t$passcount\n";
  }

}
close OUTSUM;
our $endtime=time();
#say $starttime;
#say $endtime;
printf "Done! $inputcount input line(s). %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;