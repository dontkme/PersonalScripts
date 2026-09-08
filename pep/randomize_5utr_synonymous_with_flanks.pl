#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Math::BigInt;

my $n           = 23_000;
my $id          = 'candidate';
my $input       = '';
my $seed        = undef;
my $wrap        = 0;
my $five_length = 21;
my $tail_length = 57;
my $five_flank;
my $three_flank;
my $five_flank_file;
my $three_flank_file;
my $help        = 0;

GetOptions(
    'n=i'           => \$n,
    'id=s'          => \$id,
    'input=s'       => \$input,
    'seed=i'        => \$seed,
    'wrap=i'        => \$wrap,
    'five-length=i' => \$five_length,
    'tail-length=i' => \$tail_length,
    'five-flank=s'       => \$five_flank,
    'three-flank=s'      => \$three_flank,
    'five-flank-file=s'  => \$five_flank_file,
    'three-flank-file=s' => \$three_flank_file,
    'help|h'        => \$help,
) or usage(1);

usage(0) if $help;
usage(1, '--n must be positive') if $n <= 0;
usage(1, '--wrap must be zero or positive') if $wrap < 0;
usage(1, '--five-length must be positive') if $five_length <= 0;
usage(1, '--tail-length must be a positive multiple of 3')
    if $tail_length <= 0 || $tail_length % 3;
usage(1, 'Use either --input or a command-line nucleotide sequence, not both')
    if $input ne '' && @ARGV;
usage(1, 'Use either --five-flank or --five-flank-file, not both')
    if defined $five_flank && defined $five_flank_file;
usage(1, 'Use either --three-flank or --three-flank-file, not both')
    if defined $three_flank && defined $three_flank_file;

my %aa_to_codons = (
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

my %codon_to_aa;
for my $aa (keys %aa_to_codons) {
    $codon_to_aa{$_} = $aa for @{ $aa_to_codons{$aa} };
}
my @bases = qw(A C G T);

$five_flank = defined $five_flank_file
    ? read_flank_file($five_flank_file, "5' flank")
    : normalize_flank($five_flank // '', "5' flank");
$three_flank = defined $three_flank_file
    ? read_flank_file($three_flank_file, "3' flank")
    : normalize_flank($three_flank // '', "3' flank");

my $original = read_single_sequence();
$original =~ s/\s+//g;
$original = uc $original;
$original =~ tr/U/T/;

my $expected_length = $five_length + 3 + $tail_length;
usage(1, "Input must be exactly $expected_length nt ($five_length nt + ATG + $tail_length nt); got "
    . length($original) . ' nt')
    if length($original) != $expected_length;
usage(1, 'Input contains characters other than A, C, G, T, or U')
    if $original =~ /[^ACGT]/;

my $start_codon   = substr($original, $five_length, 3);
my $original_tail = substr($original, $five_length + 3, $tail_length);
usage(1, 'Expected ATG at positions ' . ($five_length + 1) . '-'
    . ($five_length + 3) . "; found $start_codon")
    if $start_codon ne 'ATG';

my @tail_aa;
for (my $i = 0; $i < length($original_tail); $i += 3) {
    my $codon = substr($original_tail, $i, 3);
    my $aa = $codon_to_aa{$codon};
    usage(1, "Invalid codon '$codon' in the downstream segment")
        if !defined $aa;
    usage(1, 'The downstream segment contains a stop codon at relative position '
        . ($i + 1))
        if $aa eq '*';
    push @tail_aa, $aa;
}

my $tail_space = Math::BigInt->bone();
$tail_space->bmul(scalar @{ $aa_to_codons{$_} }) for @tail_aa;
usage(1, 'The downstream segment has no synonymous alternative; it contains only Met and/or Trp codons')
    if $tail_space->is_one();

# Count nucleotide strings that do not contain ATG at any position.
my $five_space = count_atg_free_sequences($five_length);
my $total_space = $five_space->copy();
$total_space->bmul($tail_space->copy()->bdec());
usage(1, "Cannot emit $n unique records; only $total_space valid sequences are possible")
    if Math::BigInt->new($n)->bcmp($total_space) > 0;

srand($seed) if defined $seed;

my %seen;
my $emitted = 0;
my $attempts = 0;
my $attempt_limit = $n * 1000;
$attempt_limit = 10_000 if $attempt_limit < 10_000;

warn "Generating $n unique records; random noncoding 5' segment=$five_length nt without ATG, "
    . "fixed start=ATG, "
    . "synonymous downstream segment=$tail_length nt.\n";
warn "Appending fixed flanks: 5'=" . length($five_flank) . " nt, 3'="
    . length($three_flank) . ' nt; output length=' .
    (length($five_flank) + $expected_length + length($three_flank)) . " nt.\n";
warn 'Downstream peptide: ' . join('', @tail_aa) . "\n";

while ($emitted < $n) {
    ++$attempts;
    die "Could not find $n unique records after $attempt_limit attempts. Try a smaller --n.\n"
        if $attempts > $attempt_limit;

    my $five = random_five_prime();
    my $tail = random_synonymous_tail();
    next if $tail eq $original_tail;

    my $sequence = $five . 'ATG' . $tail;
    next if $sequence eq $original;
    next if $seen{$sequence}++;

    ++$emitted;
    print '>', $id, '_', $emitted, "\n";
    print_wrapped($five_flank . $sequence . $three_flank);
}

exit 0;

sub random_five_prime {
    while (1) {
        my $sequence = join('', map { $bases[ int(rand(@bases)) ] } 1 .. $five_length);
        return $sequence if index($sequence, 'ATG') < 0;
    }
}

sub count_atg_free_sequences {
    my ($length) = @_;

    # States track whether the current suffix is '', 'A', or 'AT'.
    my @counts = (Math::BigInt->bone(), Math::BigInt->bzero(), Math::BigInt->bzero());
    for (1 .. $length) {
        my @next = (Math::BigInt->bzero(), Math::BigInt->bzero(), Math::BigInt->bzero());
        for my $state (0 .. 2) {
            for my $base (@bases) {
                next if $state == 2 && $base eq 'G';
                my $new_state = $base eq 'A' ? 1
                    : ($state == 1 && $base eq 'T' ? 2 : 0);
                $next[$new_state]->badd($counts[$state]);
            }
        }
        @counts = @next;
    }

    my $total = Math::BigInt->bzero();
    $total->badd($_) for @counts;
    return $total;
}

sub random_synonymous_tail {
    my $sequence = '';
    for my $aa (@tail_aa) {
        my $choices = $aa_to_codons{$aa};
        $sequence .= $choices->[ int(rand(@$choices)) ];
    }
    return $sequence;
}

sub read_single_sequence {
    return join('', @ARGV) if @ARGV;

    my $fh;
    if ($input ne '') {
        open $fh, '<', $input or die "Cannot open $input: $!\n";
    } else {
        $fh = *STDIN;
    }

    my $sequence = '';
    my $headers = 0;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            ++$headers;
            usage(1, 'Input FASTA must contain exactly one record') if $headers > 1;
            next;
        }
        $sequence .= $line;
    }
    close $fh if $input ne '';
    usage(1, 'No nucleotide sequence was provided') if $sequence =~ /^\s*$/;
    return $sequence;
}

sub read_flank_file {
    my ($path, $label) = @_;
    open my $fh, '<', $path or die "Cannot open $path: $!\n";

    my $sequence = '';
    my $headers = 0;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            ++$headers;
            usage(1, "$label file must contain at most one FASTA record")
                if $headers > 1;
            next;
        }
        $sequence .= $line;
    }
    close $fh;
    usage(1, "$label file '$path' contains no sequence")
        if $sequence =~ /^\s*$/;
    return normalize_flank($sequence, $label);
}

sub normalize_flank {
    my ($sequence, $label) = @_;
    $sequence =~ s/\s+//g;
    $sequence = uc $sequence;
    $sequence =~ tr/U/T/;
    usage(1, "$label contains characters other than A, C, G, T, or U")
        if $sequence =~ /[^ACGT]/;
    return $sequence;
}

sub print_wrapped {
    my ($sequence) = @_;
    if ($wrap == 0) {
        print $sequence, "\n";
        return;
    }
    for (my $i = 0; $i < length($sequence); $i += $wrap) {
        print substr($sequence, $i, $wrap), "\n";
    }
}

sub usage {
    my ($exit, $message) = @_;
    print STDERR "$message\n\n" if defined $message;
    print STDERR <<'USAGE';
Usage:
  perl randomize_5utr_synonymous_with_flanks.pl --input core.fa \
    --five-flank ACGT --three-flank TTAA > candidates.fa

  perl randomize_5utr_synonymous_with_flanks.pl --input core.fa \
    --five-flank-file five.txt --three-flank-file three.txt > candidates.fa

Input structure (default):
  21-nt 5' segment + ATG + 57-nt downstream coding segment = 81 nt

Options:
  --n INT             Number of unique sequences [23000]
  --id STR            FASTA record prefix [candidate]
  --input FILE        Read one raw or FASTA nucleotide sequence from FILE
  --seed INT          Make random generation reproducible
  --wrap INT          Wrap FASTA sequences at INT bases; 0 means one line [0]
  --five-length INT   Length of random noncoding 5' segment [21]
  --tail-length INT   Length of synonymous downstream segment [57]
  --five-flank SEQ    Add SEQ before every generated core sequence
  --three-flank SEQ   Add SEQ after every generated core sequence
  --five-flank-file FILE
                      Read the fixed 5' flank from a text or one-record FASTA file
  --three-flank-file FILE
                      Read the fixed 3' flank from a text or one-record FASTA file
  --help              Show this help

Rules:
  * Every base in the 5' segment is sampled independently from A/C/G/T.
    Sequences containing ATG anywhere within this segment are rejected.
  * The central ATG is preserved exactly.
  * The downstream peptide is preserved by synonymous codon sampling, and the
    generated downstream nucleotide segment must differ from the input.
  * All complete output sequences are unique and the input sequence is excluded.
  * Fixed flanks are appended after core generation and are not randomized.
  * DNA or RNA input is accepted; U is converted to T and DNA FASTA is emitted.
USAGE
    exit $exit;
}
