#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(shuffle);

my $n        = 0;
my $fraction = 0;
my $seed     = undef;
my $help     = 0;

GetOptions(
    'n=i'        => \$n,
    'fraction=f' => \$fraction,
    'seed=i'     => \$seed,
    'help|h'     => \$help,
) or usage(1);

usage(0) if $help;
usage(1, "Use exactly one of --n or --fraction") if ($n > 0) == ($fraction > 0);
usage(1, "--fraction must be between 0 and 1") if $fraction && ($fraction <= 0 || $fraction > 1);

srand($seed) if defined $seed;

if ($n > 0) {
    reservoir_sample($n);
} else {
    fraction_sample($fraction);
}

exit 0;

sub reservoir_sample {
    my ($sample_size) = @_;
    my @reservoir;
    my $seen = 0;

    read_fasta_record(sub {
        my ($record) = @_;
        ++$seen;
        if (@reservoir < $sample_size) {
            push @reservoir, $record;
            return;
        }
        my $j = int(rand($seen));
        $reservoir[$j] = $record if $j < $sample_size;
    });

    print for @reservoir;
}

sub fraction_sample {
    my ($p) = @_;
    read_fasta_record(sub {
        my ($record) = @_;
        print $record if rand() < $p;
    });
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
    $callback->($$header_ref . join('', @$seq_ref));
}

sub usage {
    my ($exit, $msg) = @_;
    print STDERR "$msg\n\n" if defined $msg;
    print STDERR <<"USAGE";
Usage:
  perl sample_fasta.pl --n 100 --seed 1 < all_nt.fa > sampled.fa
  perl sample_fasta.pl --fraction 0.01 --seed 1 < all_nt.fa > sampled.fa

Options:
  --n INT          Uniformly sample exactly INT records without replacement
  --fraction FLOAT Keep each record independently with probability FLOAT
  --seed INT       Make random sampling reproducible
  --help           Show this help
USAGE
    exit $exit;
}
