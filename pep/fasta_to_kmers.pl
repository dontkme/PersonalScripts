#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $k      = 50;
my $step   = 1;
my $wrap   = 80;
my $skip_n = 1;
my $help   = 0;

GetOptions(
    'k=i'      => \$k,
    'step=i'   => \$step,
    'wrap=i'   => \$wrap,
    'skip-n!'  => \$skip_n,
    'help|h'   => \$help,
) or usage(1);

usage(0) if $help;
usage(1, "--k must be positive") if $k <= 0;
usage(1, "--step must be positive") if $step <= 0;

my $records = 0;
my $kmers = 0;
my $skipped = 0;

read_fasta_record(sub {
    my ($header, $seq) = @_;
    ++$records;
    my $id = fasta_id($header);
    $seq = clean_seq($seq);
    return if length($seq) < $k;

    for (my $i = 0; $i <= length($seq) - $k; $i += $step) {
        my $kmer = substr($seq, $i, $k);
        if ($skip_n && $kmer =~ /[^ACGT]/) {
            ++$skipped;
            next;
        }

        my $start = $i + 1;
        my $end = $i + $k;
        print '>', $id, '|kmer|', $start, '|', $end, "\n";
        print_wrapped($kmer, $wrap);
        ++$kmers;
    }
});

my $summary = "Read $records records; emitted $kmers ${k}-mers";
$summary .= "; skipped $skipped windows containing non-ACGT bases" if $skipped;
warn "$summary.\n";
exit 0;

sub read_fasta_record {
    my ($callback) = @_;
    my $header;
    my @seq;

    while (my $line = <STDIN>) {
        if ($line =~ /^>/) {
            flush_record($callback, \$header, \@seq);
            $header = $line;
            @seq = ();
        } elsif (defined $header) {
            push @seq, $line;
        }
    }
    flush_record($callback, \$header, \@seq);
}

sub flush_record {
    my ($callback, $header_ref, $seq_ref) = @_;
    return if !defined $$header_ref;
    $callback->($$header_ref, join('', @$seq_ref));
}

sub fasta_id {
    my ($header) = @_;
    chomp $header;
    $header =~ s/^>//;
    my ($id) = split /\s+/, $header;
    usage(1, "FASTA record has an empty ID") if !defined $id || $id eq '';
    usage(1, "FASTA ID '$id' contains '|', which is reserved for k-mer metadata")
        if $id =~ /\|/;
    return $id;
}

sub clean_seq {
    my ($seq) = @_;
    $seq =~ s/\s+//g;
    return uc $seq;
}

sub print_wrapped {
    my ($seq, $wrap) = @_;
    if ($wrap <= 0) {
        print $seq, "\n";
        return;
    }
    for (my $i = 0; $i < length($seq); $i += $wrap) {
        print substr($seq, $i, $wrap), "\n";
    }
}

sub usage {
    my ($exit, $msg) = @_;
    print STDERR "$msg\n\n" if defined $msg;
    print STDERR <<"USAGE";
Usage:
  perl fasta_to_kmers.pl --k 50 --wrap 0 < candidates.fa > candidates.50mers.fa

Options:
  --k INT       K-mer length [50]
  --step INT    Step between adjacent windows [1]
  --wrap INT    Wrap output sequence lines at INT bases; 0 means no wrapping [80]
  --skip-n      Skip windows containing non-ACGT bases [default]
  --no-skip-n   Emit windows containing non-ACGT bases
  --help        Show this help

Output FASTA IDs:
  parent_id|kmer|START|END

Notes:
  START and END are 1-based coordinates in the original parent candidate.
USAGE
    exit $exit;
}
