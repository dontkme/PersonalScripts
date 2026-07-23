#!/usr/bin/env python3
"""
Compute PhyloCSF for 60nt upstream and downstream flanking regions
of MHC binding positions, using per-codon reading frame from ORF CDS.

Direction is CDS-oriented (5'→3'), accounting for strand and introns.
"""
import sys
import gzip
import argparse
import pyBigWig

# ---- ORF loading (same as add_phylocsf_to_bed.py) ----

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
            orfs[name] = {
                'chrom': cols[0], 'strand': cols[5],
                'chromStart': int(cols[1]),
                'blockSizes': [int(x) for x in cols[10].rstrip(',').split(',') if x],
                'blockStarts': [int(x) for x in cols[11].rstrip(',').split(',') if x],
            }
    return orfs

# ---- CDS ↔ genomic mapping ----

def genomic_to_cds(orf, gs, ge):
    """Map genomic interval [gs, ge) to CDS range (cds_start, cds_end).
    Returns (cds_start, cds_end) or None."""
    strand = orf['strand']
    cs = orf['chromStart']
    sizes = orf['blockSizes']
    starts = orf['blockStarts']

    if strand == '+':
        cum = 0
        for i in range(len(sizes)):
            bg = cs + starts[i]
            be = bg + sizes[i]
            if be <= gs:
                cum += sizes[i]; continue
            if bg >= ge:
                break
            if bg <= gs < be:
                cds_start = cum + (gs - bg)
            else:
                cds_start = cum
            if bg < ge <= be:
                cds_end = cum + (ge - bg)
                return (cds_start, cds_end)
            cum += sizes[i]
        return None
    else:
        cum = 0
        for i in range(len(sizes) - 1, -1, -1):
            bg = cs + starts[i]
            be = bg + sizes[i]
            if be <= gs:
                break
            if bg >= ge:
                cum += sizes[i]; continue
            # gs is leftmost genomic, ge is rightmost (exclusive)
            # On - strand: leftmost → 3' end → larger CDS offset
            #              rightmost → 5' end → smaller CDS offset
            if bg <= gs < be:
                cds_3prime = cum + (be - 1 - gs)
            else:
                cds_3prime = cum
            if bg < ge <= be:
                cds_5prime = cum + (be - 1 - (ge - 1))
                # Return (5', 3') = (smaller, larger)
                return (cds_5prime, cds_3prime + 1)  # +1 to make half-open [5', 3'+1)
            cum += sizes[i]
        return None

def cds_to_genomic_intervals(orf, cds_start, cds_end):
    """Map CDS range [cds_start, cds_end) to list of genomic intervals.
    Returns [(g_start, g_end), ...] in genomic order (left to right)."""
    strand = orf['strand']
    cs = orf['chromStart']
    sizes = orf['blockSizes']
    starts = orf['blockStarts']
    intervals = []

    if strand == '+':
        cum = 0
        for i in range(len(sizes)):
            bg = cs + starts[i]
            be = bg + sizes[i]
            block_cds_end = cum + sizes[i]
            if block_cds_end <= cds_start:
                cum += sizes[i]; continue
            if cum >= cds_end:
                break
            # Overlap
            o_cds_start = max(cds_start, cum)
            o_cds_end = min(cds_end, block_cds_end)
            o_g_start = bg + (o_cds_start - cum)
            o_g_end = bg + (o_cds_end - cum)
            intervals.append((o_g_start, o_g_end))
            cum += sizes[i]
    else:
        cum = 0
        for i in range(len(sizes) - 1, -1, -1):
            bg = cs + starts[i]
            be = bg + sizes[i]
            block_cds_end = cum + sizes[i]
            if block_cds_end <= cds_start:
                cum += sizes[i]; continue
            if cum >= cds_end:
                break
            o_cds_start = max(cds_start, cum)
            o_cds_end = min(cds_end, block_cds_end)
            # Map CDS offset to genomic (right to left)
            o_g_right = be - 1 - (o_cds_start - cum)
            o_g_left = be - 1 - (o_cds_end - 1 - cum)
            intervals.append((o_g_left, o_g_right + 1))
            cum += sizes[i]

    # Return sorted left to right
    intervals.sort()
    return intervals

# ---- PhyloCSF computation (same as add_phylocsf_to_bed.py) ----

def compute_phylocsf(orf, chrom, intervals, bw_forward, bw_reverse):
    """Compute avg PhyloCSF over a list of genomic intervals."""
    strand = orf['strand']
    cs = orf['chromStart']
    sizes = orf['blockSizes']
    starts = orf['blockStarts']
    bws = bw_forward if strand == '+' else bw_reverse
    rev_map = {0: 1, 1: 0, 2: 2}

    total = 0.0
    total_n = 0

    for (gs, ge) in intervals:
        # Which exon block?
        if strand == '+':
            cum = 0
            for i in range(len(sizes)):
                bg = cs + starts[i]
                be = bg + sizes[i]
                if be <= gs:
                    cum += sizes[i]; continue
                if bg >= ge:
                    break

                # Overlap within this block
                o_gs = max(gs, bg)
                o_ge = min(ge, be)
                all_vals = []
                for bw in bws:
                    try:
                        all_vals.append(bw.values(chrom, o_gs, o_ge))
                    except:
                        all_vals.append([float('nan')] * (o_ge - o_gs))

                n_pos = o_ge - o_gs
                for j in range(n_pos):
                    cds_offset = cum + (o_gs - bg) + j
                    codon_start_cds = (cds_offset // 3) * 3
                    codon_g_pos = bg + (codon_start_cds - cum)
                    frame_idx = codon_g_pos % 3
                    v = all_vals[frame_idx][j]
                    if v is not None and v == v:
                        total += v
                        total_n += 1

                cum += sizes[i]
        else:
            cum = 0
            for i in range(len(sizes) - 1, -1, -1):
                bg = cs + starts[i]
                be = bg + sizes[i]
                if be <= gs:
                    break
                if bg >= ge:
                    cum += sizes[i]; continue

                o_gs = max(gs, bg)
                o_ge = min(ge, be)
                all_vals = []
                for bw in bws:
                    try:
                        all_vals.append(bw.values(chrom, o_gs, o_ge))
                    except:
                        all_vals.append([float('nan')] * (o_ge - o_gs))

                n_pos = o_ge - o_gs
                for j in range(n_pos):
                    nt_g = o_gs + j
                    cds_offset = cum + (be - nt_g - 1)
                    codon_start_cds = (cds_offset // 3) * 3
                    codon_leftmost = bg + (sizes[i] - 1 - (codon_start_cds - cum + 2))
                    frame_idx = rev_map.get(codon_leftmost % 3, 0)
                    v = all_vals[frame_idx][j]
                    if v is not None and v == v:
                        total += v
                        total_n += 1

                cum += sizes[i]

    return total / total_n if total_n > 0 else 0.0


def format_intervals(intervals):
    """Format genomic intervals as string."""
    return ';'.join(f"{gs}-{ge}" for gs, ge in intervals)


def main():
    parser = argparse.ArgumentParser(
        description='PhyloCSF of flanking regions around MHC binding sites')
    parser.add_argument('bed', help='Shared BED file (.bed.gz)')
    parser.add_argument('orf_bed', help='tog.ORFs.bed')
    parser.add_argument('-o', '--output', default='flanking_phylocsf.tsv')
    parser.add_argument('--upstream', type=int, default=60)
    parser.add_argument('--downstream', type=int, default=60)
    args = parser.parse_args()

    UP = args.upstream
    DOWN = args.downstream

    print("Loading ORFs...", file=sys.stderr)
    orfs = load_orf_frames(args.orf_bed)
    print(f"  {len(orfs)} ORFs", file=sys.stderr)

    print("Loading BigWigs...", file=sys.stderr)
    bw_fwd = [pyBigWig.open(f"PhyloCSF+{i}.bw") for i in range(1, 4)]
    bw_rev = [pyBigWig.open(f"PhyloCSF-{i}.bw") for i in range(1, 4)]

    header = ['pep_name', 'orf_type', 'binding_pos', 'strand',
              'h_allele', 'h_rank', 'best_rank', 'shared_class',
              'cds_total', 'pep_cds_start', 'pep_cds_end',
              'up_intervals', 'up_phylocsf', 'up_nt',
              'down_intervals', 'down_phylocsf', 'down_nt',
              'binding_phylocsf']

    count = 0
    missing = 0
    no_flank = 0

    opener = gzip.open if args.bed.endswith('.gz') else open
    with opener(args.bed, 'rt') as fin, open(args.output, 'w') as fout:
        fout.write('\t'.join(header) + '\n')

        for line in fin:
            if line.startswith('#') or line.startswith('track'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12:
                continue

            chrom = cols[0]
            gs = int(cols[1])
            ge = int(cols[2])
            name = cols[3]
            rgb = cols[8]

            # Parse name
            parts = name.split('|')
            if len(parts) < 14:
                continue
            pep_name = '|'.join(parts[:5])
            orf_type = parts[3]
            h_allele = parts[5].replace('HUMAN:', '') if parts[5].startswith('HUMAN:') else parts[5]
            h_rank = parts[8] if 8 < len(parts) else ''

            # Shared class
            sc_map = {'255,0,0': 'red', '255,165,0': 'orange', '0,128,0': 'green', '0,0,255': 'blue'}
            shared_class = sc_map.get(rgb, rgb)

            # Best rank
            ranks = []
            for r in h_rank.split(','):
                try: ranks.append(float(r))
                except: pass
            # Mouse ranks
            for i, p in enumerate(parts):
                if p.startswith('MOUSE:') and i + 3 < len(parts):
                    for r in parts[i+3].split(','):
                        try: ranks.append(float(r))
                        except: pass
                    break
            best_rank = min(ranks) if ranks else 'NA'

            orf = orfs.get(pep_name)
            if orf is None:
                missing += 1
                continue

            # Map peptide to CDS
            cds_info = genomic_to_cds(orf, gs, ge)
            if cds_info is None:
                continue
            pep_cds_start, pep_cds_end = cds_info
            total_cds = sum(orf['blockSizes'])

            # Upstream: [max(0, pep_cds_start - UP), pep_cds_start)
            up_start = max(0, pep_cds_start - UP)
            up_end = pep_cds_start
            up_intervals = cds_to_genomic_intervals(orf, up_start, up_end) if up_start < up_end else []
            up_score = compute_phylocsf(orf, chrom, up_intervals, bw_fwd, bw_rev) if up_intervals else float('nan')
            up_nt = sum(ge2 - gs2 for gs2, ge2 in up_intervals)

            # Downstream: [pep_cds_end, min(total_cds, pep_cds_end + DOWN))
            down_start = pep_cds_end
            down_end = min(total_cds, pep_cds_end + DOWN)
            down_intervals = cds_to_genomic_intervals(orf, down_start, down_end) if down_start < down_end else []
            down_score = compute_phylocsf(orf, chrom, down_intervals, bw_fwd, bw_rev) if down_intervals else float('nan')
            down_nt = sum(ge2 - gs2 for gs2, ge2 in down_intervals)

            if up_nt == 0 and down_nt == 0:
                no_flank += 1

            # Binding site PhyloCSF
            binding_score = compute_phylocsf(orf, chrom, [(gs, ge)], bw_fwd, bw_rev)

            binding_pos = f"{chrom}:{gs}-{ge}"
            fout.write('\t'.join(map(str, [
                pep_name, orf_type, binding_pos, orf['strand'],
                h_allele, h_rank, best_rank, shared_class,
                total_cds, pep_cds_start, pep_cds_end,
                format_intervals(up_intervals), f"{up_score:.6f}", up_nt,
                format_intervals(down_intervals), f"{down_score:.6f}", down_nt,
                f"{binding_score:.6f}",
            ])) + '\n')
            count += 1
            if count % 500000 == 0:
                print(f"  {count} entries...", file=sys.stderr)

    for bw in bw_fwd + bw_rev:
        bw.close()

    print(f"Done: {count} rows, {missing} missing ORFs, {no_flank} no flank",
          file=sys.stderr)

if __name__ == '__main__':
    main()
