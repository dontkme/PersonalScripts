#!/usr/bin/env python3
"""
Add reading-frame-aware PhyloCSF scores to MHC binding BED files.

Reuses CDS frame computation from phylocsf_binding.py.
Handles both individual merged BEDs and shared human-mouse BEDs.
"""
import sys
import argparse
import pyBigWig

def load_orf_frames(bed_path):
    orfs = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12:
                continue
            name = cols[3]
            chrom = cols[0]
            strand = cols[5]
            chrom_start = int(cols[1])
            sizes = [int(x) for x in cols[10].rstrip(',').split(',') if x]
            starts = [int(x) for x in cols[11].rstrip(',').split(',') if x]
            orfs[name] = {
                'chrom': chrom, 'strand': strand,
                'chromStart': chrom_start,
                'blockSizes': sizes, 'blockStarts': starts,
            }
    return orfs

def pep_name_from_annot(name_field):
    parts = name_field.split('|')
    if len(parts) < 5:
        return name_field
    return '|'.join(parts[:5])

def compute_entry_phylocsf(orf, chrom, pep_gs, pep_ge, bw_forward, bw_reverse):
    """Compute avg PhyloCSF using per-codon frame, correct across introns.

    For each nucleotide position, determines the codon's genomic start from
    CDS offset (accounting for introns), selects the correct frame BigWig,
    and accumulates per-base scores.
    """
    strand = orf['strand']
    cs = orf['chromStart']
    sizes = orf['blockSizes']
    starts = orf['blockStarts']
    bws = bw_forward if strand == '+' else bw_reverse
    # - strand: codon_leftmost % 3 -> BigWig index
    rev_map = {0: 1, 1: 0, 2: 2}  # 0->-2, 1->-1, 2->-3

    total = 0.0
    total_n = 0

    if strand == '+':
        cum = 0
        for i in range(len(sizes)):
            size = sizes[i]
            block_g_start = cs + starts[i]
            block_g_end = block_g_start + size

            if block_g_end <= pep_gs:
                cum += size; continue
            if block_g_start >= pep_ge:
                break

            g_start = max(pep_gs, block_g_start)
            g_end = min(pep_ge, block_g_end)

            # Query all 3 BigWigs for this region
            all_vals = []
            for bw in bws:
                try:
                    all_vals.append(bw.values(chrom, g_start, g_end))
                except:
                    all_vals.append([float('nan')] * (g_end - g_start))

            n_pos = g_end - g_start
            for j in range(n_pos):
                # CDS offset of this nucleotide
                cds_offset = cum + (g_start - block_g_start) + j
                # Codon genomic start
                codon_start_cds = (cds_offset // 3) * 3
                codon_g_pos = block_g_start + (codon_start_cds - cum)
                frame_idx = codon_g_pos % 3
                v = all_vals[frame_idx][j]
                if v is not None and v == v:
                    total += v
                    total_n += 1

            cum += size
    else:
        cum = 0
        for i in range(len(sizes) - 1, -1, -1):
            size = sizes[i]
            block_g_start = cs + starts[i]
            block_g_end = block_g_start + size

            if block_g_end <= pep_gs:
                break
            if block_g_start >= pep_ge:
                cum += size; continue

            g_start = max(pep_gs, block_g_start)
            g_end = min(pep_ge, block_g_end)

            all_vals = []
            for bw in bws:
                try:
                    all_vals.append(bw.values(chrom, g_start, g_end))
                except:
                    all_vals.append([float('nan')] * (g_end - g_start))

            n_pos = g_end - g_start
            for j in range(n_pos):
                # Genomic position of this nucleotide
                nt_g = g_start + j
                # CDS offset from right end
                cds_offset = cum + (block_g_end - nt_g - 1)
                # Codon start CDS offset
                codon_start_cds = (cds_offset // 3) * 3
                # Codon's leftmost genomic position
                codon_leftmost = block_g_start + \
                    (size - 1 - (codon_start_cds - cum + 2))
                frame_idx = rev_map.get(codon_leftmost % 3, 0)
                v = all_vals[frame_idx][j]
                if v is not None and v == v:
                    total += v
                    total_n += 1

            cum += size

    return total / total_n if total_n > 0 else 0.0

def main():
    parser = argparse.ArgumentParser(
        description='Add PhyloCSF scores to MHC binding BED files')
    parser.add_argument('bed', help='Input BED file (gzipped or plain)')
    parser.add_argument('orf_bed', help='tog.ORFs.bed')
    parser.add_argument('-o', '--output', required=True, help='Output BED file (uncompressed)')
    args = parser.parse_args()

    print("Loading ORF CDS frames...", file=sys.stderr)
    orfs = load_orf_frames(args.orf_bed)
    print(f"  {len(orfs)} ORFs", file=sys.stderr)

    print("Loading PhyloCSF BigWigs...", file=sys.stderr)
    bw_fwd = [pyBigWig.open(f"PhyloCSF+{i}.bw") for i in range(1, 4)]
    bw_rev = [pyBigWig.open(f"PhyloCSF-{i}.bw") for i in range(1, 4)]

    import gzip
    opener = gzip.open if args.bed.endswith('.gz') else open
    mode = 'rt' if args.bed.endswith('.gz') else 'r'

    count = 0
    missing = 0
    with opener(args.bed, mode) as fin, open(args.output, 'w') as fout:
        for line in fin:
            if line.startswith('#') or line.startswith('track'):
                fout.write(line)
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12:
                fout.write(line + '\n')
                continue

            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])

            pep_name = pep_name_from_annot(cols[3])
            orf = orfs.get(pep_name)
            if orf is None:
                missing += 1
                cols[3] += "|phyloCSF=NA"
            else:
                score = compute_entry_phylocsf(orf, chrom, start, end, bw_fwd, bw_rev)
                cols[3] += f"|phyloCSF={score:.4f}"

            fout.write('\t'.join(cols) + '\n')
            count += 1
            if count % 1000000 == 0:
                print(f"  {count} entries...", file=sys.stderr)

    for bw in bw_fwd + bw_rev:
        bw.close()

    print(f"Done: {count} entries, {missing} missing ORFs -> {args.output}",
          file=sys.stderr)

if __name__ == '__main__':
    main()
