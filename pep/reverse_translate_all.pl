#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $id     = 'protein';
my $wrap   = 80;
my $max    = 0;
my $force  = 0;
my $help   = 0;

GetOptions(
    'id=s'    => \$id,
    'wrap=i'  => \$wrap,
    'max=i'   => \$max,
    'force!'  => \$force,
    'help|h'  => \$help,
) or usage(1);

usage(0) if $help;

my $protein = join('', @ARGV);
if ($protein eq '') {
    while (my $line = <STDIN>) {
        chomp $line;
        next if $line =~ /^>/;
        $protein .= $line;
    }
}

$protein =~ s/\s+//g;
$protein = uc $protein;
usage(1, "No protein sequence was provided") if $protein eq '';

my %codons = (
    A => [qw(GCT GCC GCA GCG)],
    R => [qw(CGT CGC CGA CGG AGA AGG)],
    N => [qw(AAT AAC)],
    D => [qw(GAT GAC)],
    C => [qw(TGT TGC)],
    Q => [qw(CAA CAG)],
    E => [qw(GAA GAG)],
    G => [qw(GGT GGC GGA GGG)],
    H => [qw(CAT CAC)],
    I => [qw(ATT ATC ATA)],
    L => [qw(TTA TTG CTT CTC CTA CTG)],
    K => [qw(AAA AAG)],
    M => [qw(ATG)],
    F => [qw(TTT TTC)],
    P => [qw(CCT CCC CCA CCG)],
    S => [qw(TCT TCC TCA TCG AGT AGC)],
    T => [qw(ACT ACC ACA ACG)],
    W => [qw(TGG)],
    Y => [qw(TAT TAC)],
    V => [qw(GTT GTC GTA GTG)],
    '*' => [qw(TAA TAG TGA)],
);

my @aa = split //, $protein;
for my $aa (@aa) {
    usage(1, "Unsupported amino acid '$aa'. Use the 20 standard one-letter amino acids; '*' is allowed for stop codons")
        if !exists $codons{$aa};
}

my $total = 1;
for my $aa (@aa) {
    $total *= scalar @{ $codons{$aa} };
    if ($max && $total > $max && !$force) {
        die "This protein expands to more than --max $max nucleotide sequences. Re-run with a larger --max or --force.\n";
    }
}

warn "Will emit $total FASTA records for protein length " . scalar(@aa) . ".\n";

my $count = 0;
emit_combinations(0, '');
exit 0;

sub emit_combinations {
    my ($pos, $seq) = @_;
    if ($pos == @aa) {
        ++$count;
        print '>', $id, '_', $count, " protein=$protein\n";
        print_wrapped($seq);
        return;
    }

    for my $codon (@{ $codons{ $aa[$pos] } }) {
        emit_combinations($pos + 1, $seq . $codon);
    }
}

sub print_wrapped {
    my ($seq) = @_;
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
  perl reverse_translate_all.pl [options] MKT
  perl reverse_translate_all.pl [options] < protein.fa > all_nt.fa

Options:
  --id STR      FASTA record prefix [protein]
  --wrap INT    Wrap sequence lines at INT bases; 0 means no wrapping [80]
  --max INT     Abort unless --force if the expansion exceeds INT records [0 = unlimited]
  --force       Allow output even when --max is exceeded
  --help        Show this help

Notes:
  Uses the standard genetic code and emits DNA codons.
  The number of records grows as the product of codon degeneracies.
USAGE
    exit $exit;
}
