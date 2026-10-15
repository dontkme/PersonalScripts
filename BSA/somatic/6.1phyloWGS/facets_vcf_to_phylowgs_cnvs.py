#!/usr/bin/env python3

import argparse
import gzip
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--facets", required=True)
parser.add_argument("--out", required=True)
parser.add_argument(
    "--purity",
    type=float,
    help="Override FACETS purity with an independently validated value"
)
args = parser.parse_args()

CNV_HEADER = (
    "chromosome\tstart\tend\tcopy_number\tminor_cn\t"
    "major_cn\tcellular_prevalence\n"
)


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def write_empty_cnvs(reason):
    """Emit a header-only CNV table so phyloWGS degrades to SNV-only.

    An empty table means "no CNV segments", which is the correct input for
    an SNV-only run -- preferable to aborting the whole step.
    """
    with open(args.out, "w") as output:
        output.write(CNV_HEADER)
    print(
        f"{reason}; wrote empty {args.out} (0 regions) so phyloWGS runs SNV-only",
        file=sys.stderr
    )


header_purity = None

with open_text(args.facets) as handle:
    for line in handle:
        if line.startswith("##purity="):
            value = line.rstrip().split("=", 1)[1]
            if value not in ("NA", ".", ""):
                header_purity = float(value)
        if line.startswith("#CHROM"):
            break

purity = args.purity if args.purity is not None else header_purity

if purity is None:
    # FACETS could not estimate purity, so no CNV segment can be scaled to a
    # cellular prevalence. Fall back to SNV-only rather than blocking the run.
    write_empty_cnvs("FACETS purity is NA")
    sys.exit(0)

if not 0 < purity <= 1:
    sys.exit("Purity must be between 0 and 1.")

if args.purity is None and abs(purity - 0.3) < 1e-10:
    sys.exit(
        "FACETS purity is exactly 0.3, which is commonly a failed estimate. "
        "Supply an independently validated --purity or run SNV-only."
    )

autosomes = {str(i) for i in range(1, 23)}
written = 0
skipped = 0

with open(args.out, "w") as output:
    output.write(
        CNV_HEADER
    )

    with open_text(args.facets) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            chrom = re.sub(r"^chr", "", fields[0], flags=re.I).upper()

            # 初次分析只使用常染色体，避免性别和单倍体校正问题
            if chrom not in autosomes:
                continue

            start = int(fields[1])
            filt = fields[6]

            info = {}
            for item in fields[7].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info[key] = value

            end_value = info.get("END")
            tcn_value = info.get("TCN_EM", ".")
            lcn_value = info.get("LCN_EM", ".")
            cf_value = info.get("CF_EM", ".")

            if end_value in (None, ".", "NA"):
                skipped += 1
                continue

            if tcn_value in (".", "NA", ""):
                skipped += 1
                continue

            end = int(end_value)
            tcn = int(float(tcn_value))

            if start >= end:
                skipped += 1
                continue

            # FACETS部分neutral片段因杂合SNP不足会缺少LCN；
            # 对明确标记为neutral且TCN=2的片段，可按1+1处理。
            if lcn_value in (".", "NA", ""):
                if tcn == 2 and filt == "neutral":
                    lcn = 1
                else:
                    skipped += 1
                    continue
            else:
                lcn = int(float(lcn_value))

            if lcn < 0 or lcn > tcn:
                skipped += 1
                continue

            major_cn = tcn - lcn

            if tcn == 2 and lcn == 1:
                # 不能使用FACETS给normal segment写的CF_EM=1；
                # 否则PhyloWGS会错误推断样本purity为100%。
                cellular_prevalence = purity
            else:
                if cf_value in (".", "NA", ""):
                    skipped += 1
                    continue

                cellular_prevalence = float(cf_value)

                # 防止小数误差或不稳定FACETS解导致CF略高于purity
                cellular_prevalence = min(cellular_prevalence, purity)

                if cellular_prevalence <= 0:
                    skipped += 1
                    continue

            output.write(
                f"{chrom}\t{start}\t{end}\t{tcn}\t{lcn}\t"
                f"{major_cn}\t{cellular_prevalence:.6f}\n"
            )
            written += 1

print(
    f"Wrote {written} regions; skipped {skipped} unusable regions",
    file=sys.stderr
)
