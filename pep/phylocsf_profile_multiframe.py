#!/usr/bin/env python3
"""
Multi-frame PhyloCSF profile: coding frame + 2 non-coding frames.

For each nucleotide position, extracts per-nt PhyloCSF from all 3 reading
frames for the correct strand. Frame0 = coding frame, Frame1/Frame2 = other
two frames (rotated from coding frame).
"""
import sys
import gzip
import argparse
import os
import pyBigWig
from flanking_phylocsf import (load_orf_frames, genomic_to_cds,
                                cds_to_genomic_intervals)


def per_nt_phylocsf_multiframe(orf, chrom, intervals, bw_forward, bw_reverse):
    """Return 3 lists of per-nt PhyloCSF scores (frame0, frame1, frame2).
    frame0 = coding frame. Order follows CDS direction (5' to 3')."""
    strand = orf['strand']
    cs = orf['chromStart']
    sizes = orf['blockSizes']
    starts = orf['blockStarts']
    bws = bw_forward if strand == '+' else bw_reverse
    rev_map = {0: 1, 1: 0, 2: 2}

    scores0, scores1, scores2 = [], [], []

    if strand == '+':
        for (gs, ge) in intervals:
            cum = 0
            for i in range(len(sizes)):
                bg = cs + starts[i]
                be = bg + sizes[i]
                if be <= gs: cum += sizes[i]; continue
                if bg >= ge: break

                o_gs = max(gs, bg)
                o_ge = min(ge, be)
                all_vals = []
                for bw in bws:
                    try:
                        all_vals.append(bw.values(chrom, o_gs, o_ge))
                    except:
                        all_vals.append([0.0] * (o_ge - o_gs))

                for j in range(o_ge - o_gs):
                    cds_offset = cum + (o_gs - bg) + j
                    codon_start_cds = (cds_offset // 3) * 3
                    codon_g_pos = bg + (codon_start_cds - cum)
                    f0 = codon_g_pos % 3
                    f1 = (f0 + 1) % 3
                    f2 = (f0 + 2) % 3
                    scores0.append(all_vals[f0][j] if (all_vals[f0][j] is not None and all_vals[f0][j] == all_vals[f0][j]) else 0.0)
                    scores1.append(all_vals[f1][j] if (all_vals[f1][j] is not None and all_vals[f1][j] == all_vals[f1][j]) else 0.0)
                    scores2.append(all_vals[f2][j] if (all_vals[f2][j] is not None and all_vals[f2][j] == all_vals[f2][j]) else 0.0)

                cum += sizes[i]
    else:
        for (gs, ge) in intervals:
            cum = 0
            for i in range(len(sizes) - 1, -1, -1):
                bg = cs + starts[i]
                be = bg + sizes[i]
                if be <= gs: break
                if bg >= ge: cum += sizes[i]; continue

                o_gs = max(gs, bg)
                o_ge = min(ge, be)
                all_vals = []
                for bw in bws:
                    try:
                        all_vals.append(bw.values(chrom, o_gs, o_ge))
                    except:
                        all_vals.append([0.0] * (o_ge - o_gs))

                n_pos = o_ge - o_gs
                for j in range(n_pos - 1, -1, -1):
                    nt_g = o_gs + j
                    cds_offset = cum + (be - nt_g - 1)
                    codon_start_cds = (cds_offset // 3) * 3
                    codon_leftmost = bg + (sizes[i] - 1 - (codon_start_cds - cum + 2))
                    f0 = rev_map.get(codon_leftmost % 3, 0)
                    f1 = (f0 + 1) % 3
                    f2 = (f0 + 2) % 3
                    scores0.append(all_vals[f0][j] if (all_vals[f0][j] is not None and all_vals[f0][j] == all_vals[f0][j]) else 0.0)
                    scores1.append(all_vals[f1][j] if (all_vals[f1][j] is not None and all_vals[f1][j] == all_vals[f1][j]) else 0.0)
                    scores2.append(all_vals[f2][j] if (all_vals[f2][j] is not None and all_vals[f2][j] == all_vals[f2][j]) else 0.0)

                cum += sizes[i]

    return scores0, scores1, scores2


def cds_range_scores_multiframe(orf, chrom, cds_start, cds_end, bw_fwd, bw_rev):
    if cds_start >= cds_end:
        return [], [], []
    intervals = cds_to_genomic_intervals(orf, cds_start, cds_end)
    if not intervals:
        return [], [], []
    return per_nt_phylocsf_multiframe(orf, chrom, intervals, bw_fwd, bw_rev)


def avg_or_na(scores):
    if not scores: return 0.0
    return sum(scores) / len(scores)


def fmt_padded(scores, pad_left, pad_right):
    parts = [''] * pad_left
    parts += [f"{v:.4f}" for v in scores]
    parts += [''] * pad_right
    return ','.join(parts)


def main():
    parser = argparse.ArgumentParser(
        description='Multi-frame PhyloCSF profile around MHC binding sites')
    parser.add_argument('bed', help='Shared BED file (.bed.gz)')
    parser.add_argument('orf_bed', help='tog.ORFs.bed')
    parser.add_argument('-o', '--output', default='phylocsf_profile_multiframe.tsv')
    parser.add_argument('--bind-len', type=int, default=27)
    parser.add_argument('--flank-up', type=int, default=60)
    parser.add_argument('--flank-down', type=int, default=60)
    args = parser.parse_args()

    BL = args.bind_len; FU = args.flank_up; FD = args.flank_down
    outdir = os.path.dirname(args.output) or '.'
    os.makedirs(outdir, exist_ok=True)

    print("Loading ORFs...", file=sys.stderr)
    orfs = load_orf_frames(args.orf_bed)
    print(f"  {len(orfs)} ORFs", file=sys.stderr)

    print("Loading BigWigs...", file=sys.stderr)
    bw_fwd = [pyBigWig.open(f"PhyloCSF+{i}.bw") for i in range(1, 4)]
    bw_rev = [pyBigWig.open(f"PhyloCSF-{i}.bw") for i in range(1, 4)]

    header = ['pep_name', 'orf_type', 'binding_pos', 'strand',
              'shared_class', 'cds_total', 'pep_cds_start', 'pep_cds_end']

    for region in ['left_bind', 'upstream', 'right_bind', 'downstream']:
        for fnum in [0, 1, 2]:
            header.append(f'{region}_frame{fnum}_scores')
        for fnum in [0, 1, 2]:
            header.append(f'{region}_frame{fnum}_avg')

    count = 0
    opener = gzip.open if args.bed.endswith('.gz') else open
    with opener(args.bed, 'rt') as fin, open(args.output, 'w') as fout:
        fout.write('\t'.join(header) + '\n')

        for line in fin:
            if line.startswith('#') or line.startswith('track'): continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12: continue

            chrom, gs, ge = cols[0], int(cols[1]), int(cols[2])
            name, rgb = cols[3], cols[8]
            parts = name.split('|')
            if len(parts) < 14: continue
            pep_name = '|'.join(parts[:5])
            orf_type = parts[3]
            sc_map = {'255,0,0': 'red', '255,165,0': 'orange',
                      '0,128,0': 'green', '0,0,255': 'blue'}
            shared_class = sc_map.get(rgb, rgb)

            orf = orfs.get(pep_name)
            if orf is None: continue

            cds_info = genomic_to_cds(orf, gs, ge)
            if cds_info is None: continue
            p_start, p_end = cds_info
            total_cds = sum(orf['blockSizes'])

            row = [pep_name, orf_type, f"{chrom}:{gs}-{ge}", orf['strand'],
                   shared_class, str(total_cds), str(p_start), str(p_end)]

            # 1. Left bind (left-aligned)
            left_end = min(p_start + BL, p_end)
            lb_s = cds_range_scores_multiframe(orf, chrom, p_start, left_end, bw_fwd, bw_rev)
            lb_pad = BL - len(lb_s[0])
            for fnum in range(3):
                row.append(fmt_padded(lb_s[fnum], 0, lb_pad))
            for fnum in range(3):
                row.append(f"{avg_or_na(lb_s[fnum]):.4f}")

            # 2. Upstream (right-aligned)
            up_s = cds_range_scores_multiframe(orf, chrom, max(0, p_start-FU), p_start, bw_fwd, bw_rev)
            up_pad = FU - len(up_s[0])
            for fnum in range(3):
                row.append(fmt_padded(up_s[fnum], up_pad, 0))
            for fnum in range(3):
                row.append(f"{avg_or_na(up_s[fnum]):.4f}")

            # 3. Right bind (right-aligned)
            right_start = max(p_start, p_end - BL)
            rb_s = cds_range_scores_multiframe(orf, chrom, right_start, p_end, bw_fwd, bw_rev)
            rb_pad = BL - len(rb_s[0])
            for fnum in range(3):
                row.append(fmt_padded(rb_s[fnum], rb_pad, 0))
            for fnum in range(3):
                row.append(f"{avg_or_na(rb_s[fnum]):.4f}")

            # 4. Downstream (left-aligned)
            down_s = cds_range_scores_multiframe(orf, chrom, p_end, min(total_cds, p_end+FD), bw_fwd, bw_rev)
            down_pad = FD - len(down_s[0])
            for fnum in range(3):
                row.append(fmt_padded(down_s[fnum], 0, down_pad))
            for fnum in range(3):
                row.append(f"{avg_or_na(down_s[fnum]):.4f}")

            fout.write('\t'.join(row) + '\n')
            count += 1
            if count % 200000 == 0:
                print(f"  {count} entries...", file=sys.stderr)

    for bw in bw_fwd + bw_rev:
        bw.close()

    print(f"Done: {count} rows -> {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
