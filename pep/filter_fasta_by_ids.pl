#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $ids_file = '';
my $wrap     = 80;
my $keep     = 0;
my $help     = 0;

GetOptions(
    'ids=s'   => \$ids_file,
    'wrap=i'  => \$wrap,
    'keep!'   => \$keep,
    'help|h'  => \$help,
) or usage(1);

usage(0) if $help;
usage(1, "--ids FILE is required") if $ids_file eq '';

my %ids;
open my $ids_fh, '<', $ids_file or die "Cannot open $ids_file: $!\n";
while (my $line = <$ids_fh>) {
    chomp $line;
    next if $line eq '';
    my ($id) = split /\t/, $line;
    $ids{$id} = 1;
}
close $ids_fh;

my $seen = 0;
my $printed = 0;
my $removed = 0;

read_fasta_record(sub {
    my ($header, $seq) = @_;
    ++$seen;
    my $id = fasta_id($header);
    my $listed = exists $ids{$id};
    my $emit = $keep ? $listed : !$listed;

    if ($emit) {
        ++$printed;
        print $header;
        print_wrapped(clean_seq($seq), $wrap);
    } else {
        ++$removed;
    }
});

warn "Read $seen records; printed $printed; removed $removed.\n";
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
  perl filter_fasta_by_ids.pl --ids reject.ids < candidates.fa > clean.fa
  perl filter_fasta_by_ids.pl --ids reject.ids --keep < candidates.fa > rejected.fa

Options:
  --ids FILE  ID list; first tab-delimited column is used
  --wrap INT  Wrap output sequence lines at INT bases; 0 means no wrapping [80]
  --keep      Keep listed IDs instead of removing them
  --help      Show this help
USAGE
    exit $exit;
}
