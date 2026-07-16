#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $k        = 50;
my $ids_file = '';
my $allow_splice = 0;
my $help     = 0;

GetOptions(
    'k=i'          => \$k,
    'ids=s'        => \$ids_file,
    'allow-splice' => \$allow_splice,
    'help|h'       => \$help,
) or usage(1);

usage(0) if $help;
usage(1, "--k must be positive") if $k <= 0;

my $ids_fh;
if ($ids_file ne '') {
    open $ids_fh, '>', $ids_file or die "Cannot write $ids_file: $!\n";
}

print join("\t", qw(parent kmer_id start end ref pos strand cigar mismatches tag)), "\n";

my %reported_parent;
my $mapped = 0;
my $exact = 0;
my $missing_mismatch_tag = 0;

while (my $line = <STDIN>) {
    next if $line =~ /^@/;
    chomp $line;
    next if $line eq '';

    my @f = split /\t/, $line;
    next if @f < 11;

    my ($qname, $flag, $rname, $pos, $mapq, $cigar) = @f[0..5];
    next if $flag & 0x4;
    next if $rname eq '*' || $cigar eq '*';
    ++$mapped;

    my ($parent, $start, $end) = parse_kmer_id($qname);
    next if $reported_parent{$parent};

    my ($mismatches, $tag_name) = mismatch_count(@f[11..$#f]);
    if (!defined $mismatches && $cigar !~ /=/) {
        ++$missing_mismatch_tag;
        next;
    }

    next if !is_full_length_exact($cigar, $mismatches, $k);

    $reported_parent{$parent} = 1;
    ++$exact;
    my $strand = ($flag & 0x10) ? '-' : '+';
    print join(
        "\t",
        $parent,
        $qname,
        defined($start) ? $start : '.',
        defined($end) ? $end : '.',
        $rname,
        $pos,
        $strand,
        $cigar,
        defined($mismatches) ? $mismatches : 0,
        defined($tag_name) ? $tag_name : 'CIGAR_EQ'
    ), "\n";
    print $ids_fh $parent, "\n" if defined $ids_fh;
}

close $ids_fh if defined $ids_fh;
warn "Scanned $mapped mapped STAR alignments; found $exact parent records with full-length exact ${k}-mer alignments.\n";
warn "Skipped $missing_mismatch_tag mapped alignments because neither NM:i nor nM:i was present and CIGAR did not use '='.\n"
    if $missing_mismatch_tag;
exit 0;

sub parse_kmer_id {
    my ($qname) = @_;
    if ($qname =~ /^(.+)\|kmer\|(\d+)\|(\d+)$/) {
        return ($1, $2, $3);
    }
    return ($qname, undef, undef);
}

sub mismatch_count {
    my (@tags) = @_;
    my $nm;
    my $nM;
    for my $tag (@tags) {
        $nm = $1 if $tag =~ /^NM:i:(\d+)$/;
        $nM = $1 if $tag =~ /^nM:i:(\d+)$/;
    }
    return ($nm, 'NM') if defined $nm;
    return ($nM, 'nM') if defined $nM;
    return (undef, undef);
}

sub is_full_length_exact {
    my ($cigar, $mismatches, $len) = @_;

    my $read_bases = 0;
    my $uses_m = 0;
    my $parsed = '';

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($op_len, $op) = ($1, $2);
        $parsed .= "$op_len$op";
        if ($op eq 'M' || $op eq '=') {
            $read_bases += $op_len;
            $uses_m = 1 if $op eq 'M';
        } elsif ($op eq 'N') {
            return 0 if !$allow_splice;
        } else {
            return 0;
        }
    }

    return 0 if $parsed ne $cigar;
    return 0 if $read_bases != $len;
    return 0 if defined($mismatches) && $mismatches != 0;
    return 0 if $uses_m && !defined($mismatches);
    return 1;
}

sub usage {
    my ($exit, $msg) = @_;
    print STDERR "$msg\n\n" if defined $msg;
    print STDERR <<"USAGE";
Usage:
  STAR ... --readFilesIn candidates.50mers.fa --outSAMtype SAM \\
    | perl star_kmer_rejects.pl --k 50 --ids reject_parent.ids > reject.tsv

  samtools view -h STAR_Aligned.out.bam \\
    | perl star_kmer_rejects.pl --k 50 --ids reject_parent.ids > reject.tsv

Options:
  --k INT     K-mer/read length that must align exactly [50]
  --ids FILE  Also write rejected parent candidate IDs, one per line
  --allow-splice
              Allow CIGAR N/skipped intron if the read consumes exactly --k
              matched bases with zero mismatches
  --help      Show this help

Exact-match rule:
  A k-mer rejects its parent only when it has a full-length, contiguous exact
  alignment: CIGAR must be 50M with NM:i:0 or nM:i:0, or 50= if the aligner
  emits '=' CIGAR. Soft clipping, insertions, deletions, skipped introns, and
  shorter alignments do not count.
  With --allow-splice, skipped introns are allowed, e.g. 25M1000N25M with
  NM:i:0 or nM:i:0 counts as a full-length exact 50 nt read match.

Recommended STAR SAM attributes:
  --outSAMattributes NH HI AS nM NM MD
USAGE
    exit $exit;
}
