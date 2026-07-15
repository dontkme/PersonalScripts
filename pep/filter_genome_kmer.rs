use std::collections::{BTreeMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::mpsc;
use std::thread;

#[derive(Clone)]
struct Config {
    genome: String,
    candidates: Option<String>,
    k: usize,
    wrap: usize,
    report: Option<String>,
    check_rc: bool,
    chunk_records: usize,
    threads: usize,
    progress: u64,
}

#[derive(Clone)]
struct Record {
    header: String,
    seq: Vec<u8>,
}

struct ChunkResult {
    index: usize,
    output: String,
    report: String,
    seen: usize,
    kept: usize,
    rejected: usize,
}

struct FastaReader<R: BufRead> {
    reader: R,
    pending_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            pending_header: None,
        }
    }

    fn read_record(&mut self) -> io::Result<Option<Record>> {
        let mut header = self.pending_header.take();
        let mut seq = Vec::new();
        let mut line = String::new();

        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }

            if line.starts_with('>') {
                let clean_header = line.trim_end_matches(&['\r', '\n'][..]).to_string();
                if header.is_some() {
                    self.pending_header = Some(clean_header);
                    break;
                }
                header = Some(clean_header);
            } else if header.is_some() {
                for b in line.bytes() {
                    if !b.is_ascii_whitespace() {
                        seq.push(b.to_ascii_uppercase());
                    }
                }
            }
        }

        Ok(header.map(|header| Record { header, seq }))
    }
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let config = parse_args()?;
    if config.genome.ends_with(".gz") {
        return Err("filter_genome_kmer.rs reads uncompressed FASTA only; use gunzip -c or the Perl version for .gz input".to_string());
    }
    if config.k == 0 || config.k > 64 {
        return Err("--k must be between 1 and 64 for u128 2-bit encoding".to_string());
    }
    if config.threads == 0 {
        return Err("--threads must be positive".to_string());
    }
    if config.threads > 1 && config.chunk_records == 0 {
        return Err("--threads requires --chunk-records INT so chunks can run in parallel".to_string());
    }

    let candidate_reader = open_candidate_reader(config.candidates.as_deref())?;
    let mut fasta = FastaReader::new(candidate_reader);
    let mut stdout = BufWriter::new(io::stdout());
    let mut report_writer = open_report(config.report.as_deref())?;

    if let Some(writer) = report_writer.as_mut() {
        writeln!(writer, "record\treason\tkmer").map_err(|e| e.to_string())?;
    }

    let (tx, rx) = mpsc::channel::<Result<ChunkResult, String>>();
    let mut active = 0usize;
    let mut chunk_index = 0usize;
    let mut next_emit = 1usize;
    let mut pending: BTreeMap<usize, ChunkResult> = BTreeMap::new();
    let mut totals = (0usize, 0usize, 0usize);

    eprintln!(
        "Running up to {} chunk{} in parallel.",
        config.threads,
        if config.threads == 1 { "" } else { "s" }
    );

    loop {
        while active >= config.threads {
            receive_one(&rx, &mut active, &mut pending)?;
            flush_ready(
                &mut pending,
                &mut next_emit,
                &mut stdout,
                report_writer.as_mut(),
                &mut totals,
            )?;
        }

        let records = read_candidate_chunk(&mut fasta, config.chunk_records)?;
        if records.is_empty() {
            break;
        }

        chunk_index += 1;
        active += 1;
        let tx = tx.clone();
        let worker_config = config.clone();
        thread::spawn(move || {
            let result = match catch_unwind(AssertUnwindSafe(|| {
                process_chunk(chunk_index, records, &worker_config)
            })) {
                Ok(result) => result,
                Err(_) => Err(format!("Chunk {chunk_index} panicked")),
            };
            let _ = tx.send(result);
        });
    }

    drop(tx);

    while active > 0 {
        receive_one(&rx, &mut active, &mut pending)?;
        flush_ready(
            &mut pending,
            &mut next_emit,
            &mut stdout,
            report_writer.as_mut(),
            &mut totals,
        )?;
    }

    stdout.flush().map_err(|e| e.to_string())?;
    if let Some(writer) = report_writer.as_mut() {
        writer.flush().map_err(|e| e.to_string())?;
    }

    eprintln!(
        "Done. Kept {} records; rejected {} of {} candidates due to shared {}-mers with the genome.",
        totals.1, totals.2, totals.0, config.k
    );
    Ok(())
}

fn receive_one(
    rx: &mpsc::Receiver<Result<ChunkResult, String>>,
    active: &mut usize,
    pending: &mut BTreeMap<usize, ChunkResult>,
) -> Result<(), String> {
    match rx.recv().map_err(|e| e.to_string())? {
        Ok(result) => {
            *active -= 1;
            pending.insert(result.index, result);
            Ok(())
        }
        Err(err) => Err(err),
    }
}

fn flush_ready(
    pending: &mut BTreeMap<usize, ChunkResult>,
    next_emit: &mut usize,
    stdout: &mut BufWriter<io::Stdout>,
    report_writer: Option<&mut BufWriter<File>>,
    totals: &mut (usize, usize, usize),
) -> Result<(), String> {
    let mut report_writer = report_writer;
    while let Some(result) = pending.remove(next_emit) {
        stdout
            .write_all(result.output.as_bytes())
            .map_err(|e| e.to_string())?;
        if let Some(writer) = report_writer.as_deref_mut() {
            writer
                .write_all(result.report.as_bytes())
                .map_err(|e| e.to_string())?;
        }
        totals.0 += result.seen;
        totals.1 += result.kept;
        totals.2 += result.rejected;
        *next_emit += 1;
    }
    Ok(())
}

fn process_chunk(
    index: usize,
    records: Vec<Record>,
    config: &Config,
) -> Result<ChunkResult, String> {
    eprintln!("Chunk {index}: loaded {} candidate records.", records.len());

    let mut query = HashSet::new();
    let mut candidate_windows = 0u64;
    for record in &records {
        for code in encoded_kmers(&record.seq, config.k) {
            candidate_windows += 1;
            query.insert(code);
            if config.check_rc {
                query.insert(revcomp_code(code, config.k));
            }
        }
    }

    eprintln!(
        "Chunk {index}: indexed {} query k-mers from {} candidate windows{}.",
        query.len(),
        candidate_windows,
        if config.check_rc {
            " including reverse complements"
        } else {
            ""
        }
    );

    let mut bad = HashSet::new();
    let (genome_windows, genome_hits) = scan_genome(config, &query, &mut bad)?;
    eprintln!(
        "Chunk {index}: scanned {genome_windows} genome windows; found {genome_hits} matching windows."
    );

    let mut output = String::new();
    let mut report = String::new();
    let mut kept = 0usize;
    let mut rejected = 0usize;

    for record in &records {
        if let Some(hit) = first_bad_kmer(&record.seq, config.k, &bad) {
            rejected += 1;
            if config.report.is_some() {
                report.push_str(record_id_for_report(&record.header));
                report.push('\t');
                report.push_str(&format!("genome_{}mer\t{}\n", config.k, decode_kmer(hit, config.k)));
            }
        } else {
            kept += 1;
            output.push_str(&record.header);
            output.push('\n');
            push_wrapped(&mut output, &record.seq, config.wrap);
        }
    }

    eprintln!("Chunk {index}: kept {kept} records; rejected {rejected} records.");

    Ok(ChunkResult {
        index,
        output,
        report,
        seen: records.len(),
        kept,
        rejected,
    })
}

fn scan_genome(
    config: &Config,
    query: &HashSet<u128>,
    bad: &mut HashSet<u128>,
) -> Result<(u64, u64), String> {
    let file = File::open(&config.genome).map_err(|e| format!("Cannot open {}: {e}", config.genome))?;
    let reader = BufReader::new(file);
    let mut windows = 0u64;
    let mut hits = 0u64;
    let mut next_progress = config.progress;
    let mask = kmer_mask(config.k);
    let mut code = 0u128;
    let mut valid = 0usize;

    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('>') {
            code = 0;
            valid = 0;
            continue;
        }

        for b in line.bytes() {
            if let Some(bits) = base_bits(b) {
                code = ((code << 2) | bits as u128) & mask;
                valid += 1;
                if valid >= config.k {
                    windows += 1;
                    if query.contains(&code) {
                        hits += 1;
                        bad.insert(code);
                        if config.check_rc {
                            bad.insert(revcomp_code(code, config.k));
                        }
                    }
                    if config.progress > 0 && windows >= next_progress {
                        eprintln!("  scanned {windows} genome windows; current hits {hits}.");
                        next_progress += config.progress;
                    }
                }
            } else {
                code = 0;
                valid = 0;
            }
        }
    }

    Ok((windows, hits))
}

fn encoded_kmers(seq: &[u8], k: usize) -> Vec<u128> {
    let mut out = Vec::new();
    let mask = kmer_mask(k);
    let mut code = 0u128;
    let mut valid = 0usize;

    for &b in seq {
        if let Some(bits) = base_bits(b) {
            code = ((code << 2) | bits as u128) & mask;
            valid += 1;
            if valid >= k {
                out.push(code);
            }
        } else {
            code = 0;
            valid = 0;
        }
    }

    out
}

fn first_bad_kmer(seq: &[u8], k: usize, bad: &HashSet<u128>) -> Option<u128> {
    let mask = kmer_mask(k);
    let mut code = 0u128;
    let mut valid = 0usize;

    for &b in seq {
        if let Some(bits) = base_bits(b) {
            code = ((code << 2) | bits as u128) & mask;
            valid += 1;
            if valid >= k && bad.contains(&code) {
                return Some(code);
            }
        } else {
            code = 0;
            valid = 0;
        }
    }

    None
}

fn base_bits(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn kmer_mask(k: usize) -> u128 {
    if k == 64 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    }
}

fn revcomp_code(mut code: u128, k: usize) -> u128 {
    let mut rc = 0u128;
    for _ in 0..k {
        let bits = code & 3;
        rc = (rc << 2) | (3 - bits);
        code >>= 2;
    }
    rc
}

fn decode_kmer(mut code: u128, k: usize) -> String {
    let mut bases = vec![b'A'; k];
    for i in (0..k).rev() {
        bases[i] = match code & 3 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        code >>= 2;
    }
    String::from_utf8(bases).expect("encoded k-mer should be valid ASCII")
}

fn push_wrapped(out: &mut String, seq: &[u8], wrap: usize) {
    if wrap == 0 {
        out.push_str(&String::from_utf8_lossy(seq));
        out.push('\n');
        return;
    }

    for chunk in seq.chunks(wrap) {
        out.push_str(&String::from_utf8_lossy(chunk));
        out.push('\n');
    }
}

fn record_id_for_report(header: &str) -> &str {
    header.strip_prefix('>').unwrap_or(header)
}

fn read_candidate_chunk<R: BufRead>(
    fasta: &mut FastaReader<R>,
    limit: usize,
) -> Result<Vec<Record>, String> {
    let mut records = Vec::new();
    while limit == 0 || records.len() < limit {
        match fasta.read_record().map_err(|e| e.to_string())? {
            Some(record) => records.push(record),
            None => break,
        }
    }
    Ok(records)
}

fn open_candidate_reader(path: Option<&str>) -> Result<Box<dyn BufRead>, String> {
    if let Some(path) = path {
        if path.ends_with(".gz") {
            return Err("candidate .gz input is not supported by the Rust version; decompress first".to_string());
        }
        let file = File::open(path).map_err(|e| format!("Cannot open {path}: {e}"))?;
        Ok(Box::new(BufReader::new(file)))
    } else {
        Ok(Box::new(BufReader::new(io::stdin())))
    }
}

fn open_report(path: Option<&str>) -> Result<Option<BufWriter<File>>, String> {
    match path {
        Some(path) => {
            let file = File::create(path).map_err(|e| format!("Cannot write {path}: {e}"))?;
            Ok(Some(BufWriter::new(file)))
        }
        None => Ok(None),
    }
}

fn parse_args() -> Result<Config, String> {
    let mut config = Config {
        genome: String::new(),
        candidates: None,
        k: 50,
        wrap: 80,
        report: None,
        check_rc: true,
        chunk_records: 0,
        threads: 1,
        progress: 100_000_000,
    };

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--genome" => config.genome = next_value(&mut args, "--genome")?,
            "--candidates" => config.candidates = Some(next_value(&mut args, "--candidates")?),
            "--k" => config.k = parse_usize(&next_value(&mut args, "--k")?, "--k")?,
            "--wrap" => config.wrap = parse_usize(&next_value(&mut args, "--wrap")?, "--wrap")?,
            "--report" => config.report = Some(next_value(&mut args, "--report")?),
            "--no-rc" => config.check_rc = false,
            "--chunk-records" => {
                config.chunk_records = parse_usize(&next_value(&mut args, "--chunk-records")?, "--chunk-records")?
            }
            "--threads" => config.threads = parse_usize(&next_value(&mut args, "--threads")?, "--threads")?,
            "--progress" => config.progress = parse_u64(&next_value(&mut args, "--progress")?, "--progress")?,
            "--help" | "-h" => {
                print_usage();
                std::process::exit(0);
            }
            _ => return Err(format!("Unknown option: {arg}\n\n{}", usage_text())),
        }
    }

    if config.genome.is_empty() {
        return Err(format!("--genome FILE is required\n\n{}", usage_text()));
    }

    Ok(config)
}

fn next_value<I>(args: &mut I, name: &str) -> Result<String, String>
where
    I: Iterator<Item = String>,
{
    args.next()
        .ok_or_else(|| format!("{name} requires a value"))
}

fn parse_usize(value: &str, name: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .map_err(|_| format!("{name} requires a non-negative integer"))
}

fn parse_u64(value: &str, name: &str) -> Result<u64, String> {
    value
        .parse::<u64>()
        .map_err(|_| format!("{name} requires a non-negative integer"))
}

fn print_usage() {
    eprint!("{}", usage_text());
}

fn usage_text() -> &'static str {
    "Usage:
  filter_genome_kmer_rs --genome hg38.fa --k 50 < candidates.fa > kept.fa
  filter_genome_kmer_rs --genome hg38.fa --candidates candidates.fa --threads 8 --chunk-records 2500 > kept.fa

Options:
  --genome FILE         Genome FASTA to scan; uncompressed FASTA only
  --candidates FILE     Candidate FASTA [STDIN]
  --k INT               Reject a candidate with any exact INT-nt genome match [50]
  --wrap INT            Wrap output sequence lines at INT bases; 0 means no wrapping [80]
  --report FILE         Write rejected record IDs and matching k-mers to TSV
  --no-rc               Do not check reverse-complement matches [default checks both strands]
  --chunk-records INT   Process INT candidates per chunk [0 = all at once]
  --threads INT         Run up to INT chunks in parallel; requires --chunk-records [1]
  --progress INT        Print progress every INT genome windows; 0 disables [100000000]
  --help                Show this help

Notes:
  Uses 2-bit rolling k-mers in u128, so --k must be <= 64.
  The genome is scanned once per candidate chunk. Larger chunks are faster but
  use more memory; smaller chunks rescan the genome more times.
"
}
