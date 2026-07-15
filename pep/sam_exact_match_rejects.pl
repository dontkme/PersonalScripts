#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $k        = 50;
my $ids_file = '';
my $help     = 0;

GetOptions(
    'k=i'     => \$k,
    'ids=s'   => \$ids_file,
    'help|h'  => \$help,
) or usage(1);

usage(0) if $help;
usage(1, "--k must be positive") if $k <= 0;

my $ids_fh;
if ($ids_file ne '') {
    open $ids_fh, '>', $ids_file or die "Cannot write $ids_file: $!\n";
}

print join("\t", qw(record ref pos strand exact_len method cigar md)), "\n";

my %reported;
my $alignments = 0;
my $rejected = 0;
my $missing_md_for_m = 0;

while (my $line = <STDIN>) {
    next if $line =~ /^@/;
    chomp $line;
    next if $line eq '';

    my @f = split /\t/, $line;
    next if @f < 11;

    my ($qname, $flag, $rname, $pos, $mapq, $cigar) = @f[0..5];
    next if $flag & 0x4;
    next if $rname eq '*' || $cigar eq '*';
    next if $reported{$qname};
    ++$alignments;

    my $md = '';
    for my $field (@f[11..$#f]) {
        if ($field =~ /^MD:Z:(.+)$/) {
            $md = $1;
            last;
        }
    }

    my ($best_len, $best_pos, $method, $needed_md) = longest_exact_match($cigar, $md, $pos);
    ++$missing_md_for_m if $needed_md && $md eq '';
    next if $best_len < $k;

    $reported{$qname} = 1;
    ++$rejected;
    my $strand = ($flag & 0x10) ? '-' : '+';
    print join("\t", $qname, $rname, $best_pos, $strand, $best_len, $method, $cigar, ($md eq '' ? '.' : $md)), "\n";
    print $ids_fh $qname, "\n" if defined $ids_fh;
}

close $ids_fh if defined $ids_fh;
warn "Scanned $alignments mapped alignments; found $rejected records with exact matches >= $k nt.\n";
warn "Saw CIGAR M operations without MD tags in $missing_md_for_m alignments; exact length inside M could not be evaluated.\n"
    if $missing_md_for_m;
exit 0;

sub longest_exact_match {
    my ($cigar, $md, $start_pos) = @_;
    my @md_status = $md eq '' ? () : parse_md($md);
    my $uses_eqx = $cigar =~ /[=X]/ ? 1 : 0;
    my $needed_md = $cigar =~ /M/ ? 1 : 0;
    my $method = $uses_eqx ? 'cigar_eqx' : ($md ne '' ? 'md' : 'unknown');

    my $ref_pos = $start_pos;
    my $run_len = 0;
    my $run_start = $start_pos;
    my $best_len = 0;
    my $best_pos = $start_pos;
    my $md_i = 0;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($len, $op) = ($1, $2);

        if ($op eq '=') {
            if ($run_len == 0) {
                $run_start = $ref_pos;
            }
            $run_len += $len;
            ($best_len, $best_pos) = ($run_len, $run_start) if $run_len > $best_len;
            $ref_pos += $len;
        } elsif ($op eq 'X') {
            $run_len = 0;
            $ref_pos += $len;
        } elsif ($op eq 'M') {
            if ($md eq '') {
                $run_len = 0;
                $ref_pos += $len;
                next;
            }
            for (my $i = 0; $i < $len; ++$i) {
                my $status = $md_status[$md_i++] // 'X';
                if ($status eq 'M') {
                    if ($run_len == 0) {
                        $run_start = $ref_pos;
                    }
                    ++$run_len;
                    ($best_len, $best_pos) = ($run_len, $run_start) if $run_len > $best_len;
                } else {
                    $run_len = 0;
                }
                ++$ref_pos;
            }
        } elsif ($op eq 'D' || $op eq 'N') {
            $run_len = 0;
            for (my $i = 0; $i < $len && $md_i < @md_status && $md_status[$md_i] eq 'D'; ++$i) {
                ++$md_i;
            }
            $ref_pos += $len;
        } elsif ($op eq 'I' || $op eq 'S' || $op eq 'H' || $op eq 'P') {
            $run_len = 0;
        }
    }

    return ($best_len, $best_pos, $method, $needed_md);
}

sub parse_md {
    my ($md) = @_;
    my @status;

    while ($md ne '') {
        if ($md =~ s/^(\d+)//) {
            push @status, ('M') x $1;
        } elsif ($md =~ s/^\^([A-Za-z]+)//) {
            push @status, ('D') x length($1);
        } elsif ($md =~ s/^([A-Za-z])//) {
            push @status, 'X';
        } else {
            die "Cannot parse MD tag near '$md'\n";
        }
    }

    return @status;
}

sub usage {
    my ($exit, $msg) = @_;
    print STDERR "$msg\n\n" if defined $msg;
    print STDERR <<"USAGE";
Usage:
  bwa mem -a hg38.fa candidates.fa | perl sam_exact_match_rejects.pl --k 50 --ids reject.ids > reject.tsv
  minimap2 -a --eqx hg38.mmi candidates.fa | perl sam_exact_match_rejects.pl --k 50 --ids reject.ids > reject.tsv

Options:
  --k INT     Report records with any exact aligned segment at least INT nt [50]
  --ids FILE  Also write rejected record IDs, one per line
  --help      Show this help

Notes:
  CIGAR '=' operations are treated as exact matches directly.
  CIGAR 'M' operations require an MD:Z tag, as produced by BWA-MEM, to separate
  matches from mismatches. Alignments without MD tags and without '=' CIGAR
  cannot be evaluated safely.
USAGE
    exit $exit;
}
