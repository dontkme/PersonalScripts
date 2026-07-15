#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempdir);
use IO::Uncompress::Gunzip qw($GunzipError);

my $genome        = '';
my $candidates    = '';
my $k             = 50;
my $wrap          = 80;
my $report        = '';
my $check_rc      = 1;
my $chunk_records = 0;
my $threads       = 1;
my $progress      = 100_000_000;
my $help          = 0;

GetOptions(
    'genome=s'        => \$genome,
    'candidates=s'    => \$candidates,
    'k=i'             => \$k,
    'wrap=i'          => \$wrap,
    'report=s'        => \$report,
    'rc!'             => \$check_rc,
    'chunk-records=i' => \$chunk_records,
    'threads=i'       => \$threads,
    'progress=i'      => \$progress,
    'help|h'          => \$help,
) or usage(1);

usage(0) if $help;
usage(1, "--genome FILE is required") if $genome eq '';
usage(1, "--k must be positive") if $k <= 0;
usage(1, "--chunk-records must be non-negative") if $chunk_records < 0;
usage(1, "--threads must be positive") if $threads <= 0;
usage(1, "--threads requires --chunk-records INT so chunks can run in parallel")
    if $threads > 1 && $chunk_records == 0;

if ($threads == 1) {
    run_sequential();
} else {
    run_parallel();
}

exit 0;

sub process_chunk {
    my ($chunk_index, $records, $out_fh, $chunk_report_fh) = @_;
    my $seen = scalar @$records;
    warn "Chunk $chunk_index: loaded $seen candidate records.\n";

    my %query_kmer;
    my $candidate_windows = build_query_kmers($records, \%query_kmer);
    warn "Chunk $chunk_index: indexed " . scalar(keys %query_kmer)
        . " query k-mers from $candidate_windows candidate windows"
        . ($check_rc ? " including reverse complements" : "") . ".\n";

    my %bad_kmer;
    my ($genome_windows, $genome_hits) = scan_genome(\%query_kmer, \%bad_kmer);
    warn "Chunk $chunk_index: scanned $genome_windows genome windows; found $genome_hits matching windows.\n";

    my ($kept, $rejected) = emit_filtered_records($records, \%bad_kmer, $out_fh, $chunk_report_fh);
    warn "Chunk $chunk_index: kept $kept records; rejected $rejected records.\n";
    return ($seen, $kept, $rejected, $candidate_windows, $genome_windows, $genome_hits);
}

sub run_sequential {
    my $report_fh = open_report($report, 1);
    my $candidate_fh = open_fasta($candidates);
    my $total_seen = 0;
    my $total_kept = 0;
    my $total_rejected = 0;
    my $chunk_index = 0;

    while (1) {
        my @records = read_candidate_chunk($candidate_fh, $chunk_records);
        last if !@records;
        ++$chunk_index;

        my ($seen, $kept, $rejected) = process_chunk($chunk_index, \@records, *STDOUT, $report_fh);
        $total_seen += $seen;
        $total_kept += $kept;
        $total_rejected += $rejected;
    }

    close $report_fh if defined $report_fh;
    warn "Done. Kept $total_kept records; rejected $total_rejected of $total_seen candidates due to shared ${k}-mers with the genome.\n";
}

sub run_parallel {
    my $tmpdir = tempdir('filter_genome_kmer.XXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $candidate_fh = open_fasta($candidates);
    my @chunks;
    my %active;
    my $chunk_index = 0;

    warn "Running up to $threads chunks in parallel; temporary files are in $tmpdir.\n";

    while (1) {
        my @records = read_candidate_chunk($candidate_fh, $chunk_records);
        last if !@records;
        ++$chunk_index;

        wait_for_one_child(\%active) while scalar(keys %active) >= $threads;

        my $chunk = {
            index  => $chunk_index,
            out    => "$tmpdir/chunk_${chunk_index}.fa",
            report => "$tmpdir/chunk_${chunk_index}.tsv",
            stats  => "$tmpdir/chunk_${chunk_index}.stats",
        };

        my $pid = fork();
        die "fork failed: $!\n" if !defined $pid;

        if ($pid == 0) {
            open my $out_fh, '>', $chunk->{out} or die "Cannot write $chunk->{out}: $!\n";
            my $chunk_report_fh = open_report($chunk->{report}, 0);
            my @stats = process_chunk($chunk_index, \@records, $out_fh, $chunk_report_fh);
            close $out_fh;
            close $chunk_report_fh if defined $chunk_report_fh;

            open my $stats_fh, '>', $chunk->{stats} or die "Cannot write $chunk->{stats}: $!\n";
            print $stats_fh join("\t", @stats), "\n";
            close $stats_fh;
            exit 0;
        }

        $chunk->{pid} = $pid;
        $active{$pid} = $chunk;
        push @chunks, $chunk;
    }

    wait_for_one_child(\%active) while %active;

    my $final_report_fh = open_report($report, 1);
    my ($total_seen, $total_kept, $total_rejected) = (0, 0, 0);

    for my $chunk (sort { $a->{index} <=> $b->{index} } @chunks) {
        my @stats = read_stats($chunk->{stats});
        $total_seen += $stats[0];
        $total_kept += $stats[1];
        $total_rejected += $stats[2];
        copy_file_to_fh($chunk->{out}, *STDOUT);
        copy_file_to_fh($chunk->{report}, $final_report_fh) if defined $final_report_fh;
    }

    close $final_report_fh if defined $final_report_fh;
    warn "Done. Kept $total_kept records; rejected $total_rejected of $total_seen candidates due to shared ${k}-mers with the genome.\n";
}

sub wait_for_one_child {
    my ($active) = @_;
    my $pid = wait();
    die "wait failed: $!\n" if $pid < 0;
    my $status = $?;
    my $chunk = delete $active->{$pid};
    die "Internal error: waited for unknown child $pid\n" if !defined $chunk;
    die "Chunk $chunk->{index} failed with status $status\n" if $status != 0;
}

sub build_query_kmers {
    my ($records, $query_kmer) = @_;
    my $windows = 0;

    for my $record (@$records) {
        my ($header, $seq) = @$record;
        my @kmers = kmers($seq, $k);
        $windows += @kmers;
        for my $kmer (@kmers) {
            $query_kmer->{$kmer} = 1;
            $query_kmer->{ revcomp($kmer) } = 1 if $check_rc;
        }
    }

    return $windows;
}

sub scan_genome {
    my ($query_kmer, $bad_kmer) = @_;
    my $fh = open_fasta($genome);
    my $carry = '';
    my $windows = 0;
    my $hits = 0;
    my $next_progress = $progress;

    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            $carry = '';
            next;
        }

        $line = uc $line;
        my $last_end = 0;
        while ($line =~ /([ACGT]+)/g) {
            $carry = '' if $-[0] > $last_end;
            my $text = $carry . $1;
            my $limit = length($text) - $k;

            for (my $i = 0; $i <= $limit; ++$i) {
                my $kmer = substr($text, $i, $k);
                ++$windows;
                if (exists $query_kmer->{$kmer}) {
                    ++$hits;
                    $bad_kmer->{$kmer} = 1;
                    $bad_kmer->{ revcomp($kmer) } = 1 if $check_rc;
                }
            }

            if (length($text) >= $k - 1) {
                $carry = substr($text, -($k - 1));
            } else {
                $carry = $text;
            }

            if ($progress > 0 && $windows >= $next_progress) {
                warn "  scanned $windows genome windows; current hits $hits.\n";
                $next_progress += $progress;
            }

            $last_end = $+[0];
        }

        $carry = '' if $last_end < length($line);
    }

    return ($windows, $hits);
}

sub emit_filtered_records {
    my ($records, $bad_kmer, $out_fh, $chunk_report_fh) = @_;
    my $kept = 0;
    my $rejected = 0;

    for my $record (@$records) {
        my ($header, $seq) = @$record;
        my $hit = first_bad_kmer($seq, $bad_kmer);

        if (defined $hit) {
            ++$rejected;
            if (defined $chunk_report_fh) {
                my $clean_header = $header;
                chomp $clean_header;
                $clean_header =~ s/^>//;
                print $chunk_report_fh join("\t", $clean_header, "genome_${k}mer", $hit), "\n";
            }
            next;
        }

        ++$kept;
        print $out_fh $header;
        print_wrapped($out_fh, $seq, $wrap);
    }

    return ($kept, $rejected);
}

sub first_bad_kmer {
    my ($seq, $bad_kmer) = @_;
    for my $kmer (kmers($seq, $k)) {
        return $kmer if exists $bad_kmer->{$kmer};
    }
    return undef;
}

sub kmers {
    my ($seq, $len) = @_;
    my @out;
    return @out if length($seq) < $len;

    for (my $i = 0; $i <= length($seq) - $len; ++$i) {
        my $kmer = substr($seq, $i, $len);
        next if $kmer =~ /[^ACGT]/;
        push @out, $kmer;
    }

    return @out;
}

sub read_candidate_chunk {
    my ($fh, $limit) = @_;
    my @records;

    while ($limit == 0 || @records < $limit) {
        my $record = read_fasta_record($fh);
        last if !defined $record;
        push @records, $record;
    }

    return @records;
}

{
    my $pending_header;

    sub read_fasta_record {
        my ($fh) = @_;
        my $header = $pending_header;
        $pending_header = undef;
        my @seq;

        while (my $line = <$fh>) {
            if ($line =~ /^>/) {
                if (defined $header) {
                    $pending_header = $line;
                    return [$header, clean_seq(join('', @seq))];
                }
                $header = $line;
            } elsif (defined $header) {
                push @seq, $line;
            }
        }

        return undef if !defined $header;
        return [$header, clean_seq(join('', @seq))];
    }
}

sub clean_seq {
    my ($seq) = @_;
    $seq =~ s/\s+//g;
    return uc $seq;
}

sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGT/TGCA/;
    return $seq;
}

sub print_wrapped {
    my ($fh, $seq, $wrap) = @_;
    if ($wrap <= 0) {
        print $fh $seq, "\n";
        return;
    }
    for (my $i = 0; $i < length($seq); $i += $wrap) {
        print $fh substr($seq, $i, $wrap), "\n";
    }
}

sub open_report {
    my ($path, $include_header) = @_;
    return undef if $path eq '';
    open my $fh, '>', $path or die "Cannot write $path: $!\n";
    print $fh join("\t", qw(record reason kmer)), "\n" if $include_header;
    return $fh;
}

sub read_stats {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot open $path: $!\n";
    my $line = <$fh>;
    close $fh;
    die "Missing stats in $path\n" if !defined $line;
    chomp $line;
    return split /\t/, $line;
}

sub copy_file_to_fh {
    my ($path, $out_fh) = @_;
    open my $in_fh, '<', $path or die "Cannot open $path: $!\n";
    while (my $line = <$in_fh>) {
        print $out_fh $line;
    }
    close $in_fh;
}

sub open_fasta {
    my ($path) = @_;
    if ($path eq '') {
        return *STDIN;
    }

    if ($path =~ /\.gz$/) {
        my $fh = IO::Uncompress::Gunzip->new($path)
            or die "Cannot open gzip file $path: $GunzipError\n";
        return $fh;
    }

    open my $fh, '<', $path or die "Cannot open $path: $!\n";
    return $fh;
}

sub usage {
    my ($exit, $msg) = @_;
    print STDERR "$msg\n\n" if defined $msg;
    print STDERR <<"USAGE";
Usage:
  perl filter_genome_kmer.pl --genome hg38.fa --k 50 < candidates.fa > kept.fa
  perl filter_genome_kmer.pl --genome hg38.fa.gz --candidates candidates.fa --report rejected.tsv > kept.fa

Options:
  --genome FILE         Genome FASTA to scan; .gz is supported
  --candidates FILE     Candidate FASTA [STDIN]
  --k INT               Reject a candidate with any exact INT-nt genome match [50]
  --wrap INT            Wrap output sequence lines at INT bases; 0 means no wrapping [80]
  --report FILE         Write rejected record IDs and matching k-mers to TSV
  --no-rc               Do not check reverse-complement matches [default checks both strands]
  --chunk-records INT   Process INT candidates per genome scan [0 = all at once]
  --threads INT         Run up to INT chunks in parallel; requires --chunk-records [1]
  --progress INT        Print progress every INT genome windows; 0 disables [100000000]
  --help                Show this help

Notes:
  This is an exact k-mer filter, not a gapped or approximate aligner. It loads
  candidate k-mers into memory, streams through the genome, then emits only
  candidates that do not share any exact k-mer with either strand of the genome.
  Parallel mode makes each worker scan the genome for one candidate chunk, then
  merges chunk outputs in input order.
USAGE
    exit $exit;
}
