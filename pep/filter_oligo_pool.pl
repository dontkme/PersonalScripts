#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $k       = 50;
my $max     = 0;
my $seed    = undef;
my $shuffle = 0;
my $wrap    = 80;
my $help    = 0;

GetOptions(
    'k=i'       => \$k,
    'max=i'     => \$max,
    'seed=i'    => \$seed,
    'shuffle!'  => \$shuffle,
    'wrap=i'    => \$wrap,
    'help|h'    => \$help,
) or usage(1);

usage(0) if $help;
usage(1, "--k must be positive") if $k <= 0;
srand($seed) if defined $seed;

my @records;
read_fasta_record(sub {
    my ($header, $seq) = @_;
    $seq =~ s/\s+//g;
    $seq = uc $seq;
    push @records, [$header, $seq];
});

@records = fisher_yates(@records) if $shuffle;

my %selected_kmer;
my $selected = 0;
my $rejected = 0;

for my $record (@records) {
    my ($header, $seq) = @$record;
    my @kmers = unique_kmers($seq, $k);

    my $conflict = 0;
    for my $kmer (@kmers) {
        if (exists $selected_kmer{$kmer}) {
            $conflict = 1;
            last;
        }
    }

    if ($conflict) {
        ++$rejected;
        next;
    }

    print $header;
    print_wrapped($seq, $wrap);
    ++$selected;
    $selected_kmer{$_} = 1 for @kmers;
    last if $max && $selected >= $max;
}

warn "Selected $selected records; rejected $rejected records due to shared ${k}-mers.\n";
exit 0;

sub unique_kmers {
    my ($seq, $len) = @_;
    my %seen;
    return ($seq) if length($seq) < $len;
    for (my $i = 0; $i <= length($seq) - $len; ++$i) {
        $seen{ substr($seq, $i, $len) } = 1;
    }
    return keys %seen;
}

sub fisher_yates {
    my (@items) = @_;
    for (my $i = @items - 1; $i > 0; --$i) {
        my $j = int(rand($i + 1));
        @items[$i, $j] = @items[$j, $i];
    }
    return @items;
}

sub read_fasta_record {
    my ($callback) = @_;
    my $header = undef;
    my @seq;

    while (my $line = <STDIN>) {
        if ($line =~ /^>/) {
            flush_record($callback, \$header, \@seq);
            $header = $line;
            @seq = ();
        } else {
            push @seq, $line if defined $header;
        }
    }
    flush_record($callback, \$header, \@seq);
}

sub flush_record {
    my ($callback, $header_ref, $seq_ref) = @_;
    return if !defined $$header_ref;
    my $seq = join('', @$seq_ref);
    $callback->($$header_ref, $seq);
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
  perl filter_oligo_pool.pl --k 50 --shuffle --seed 1 --max 1000 < sampled.fa > pool.fa

Options:
  --k INT       Reject a candidate if it shares any exact INT-nt word with a selected sequence [50]
  --max INT     Stop after selecting INT records [0 = no limit]
  --shuffle     Randomize input order before greedy selection
  --seed INT    Make shuffling reproducible
  --wrap INT    Wrap output sequence lines at INT bases; 0 means no wrapping [80]
  --help        Show this help

Strategy:
  Greedily accepts a sequence only if none of its k-mers have already appeared
  in the accepted set. This enforces "no exact matching 50 nt window" among
  selected oligos when --k 50 is used.
USAGE
    exit $exit;
}
