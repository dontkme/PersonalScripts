#!/usr/bin/env python3
"""Prepare pVACseq inputs for step 8 (VacQ).

Two independent subcommands, both run with the base python3 (the pvactools
env activation happens separately, and only around VEP/pVACseq):

  vacq_prepare.py hla   <R1_bestguess_G.txt>    -> HLA-A*11:50,HLA-A*11:01,...
  vacq_prepare.py purity <facet_out.vcf.gz>     -> "--tumor-purity 0.449" or ""

HLA: class I loci A/B/C only (NetMHCpan MHC class I cannot use class II),
restricted to alleles HLA-LA marked perfectG=1, normalised from
G-group names to the 2-field form pVACseq expects.

Purity: omits the flag entirely when FACETS' estimate is unusable -- NA,
missing, out of (0, 1], or exactly 0.3 (FACETS' failed-estimate sentinel).
Better to let pVACseq use its own default than feed it a known-bad value.
"""

import argparse
import gzip
import re
import sys

CLASS_I_LOCI = ("A", "B", "C")

# Expression-level suffixes on HLA allele names: G (G-group), Q (questionable
# expression), N (null), S (secreted), L (low expression).
EXPR_SUFFIX = re.compile(r"[A-Za-z]+$")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def format_allele(raw):
    """'A*11:01:01G' -> 'HLA-A*11:01'; returns None if unparseable."""
    allele = raw.strip()
    if not allele or "*" not in allele:
        return None
    gene, rest = allele.split("*", 1)
    fields = rest.split(":")
    if len(fields) < 2:
        return None
    fields[-1] = EXPR_SUFFIX.sub("", fields[-1])
    if not fields[0] or not fields[-1]:
        return None
    return f"HLA-{gene}*{fields[0]}:{fields[1]}"


def cmd_hla(args):
    rows = []
    with open_text(args.bestguess) as handle:
        header = None
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))

    if not rows:
        sys.exit(f"No rows parsed from {args.bestguess}")

    by_locus = {}
    for row in rows:
        locus = (row.get("Locus") or "").strip()
        if locus not in CLASS_I_LOCI:
            continue
        if (row.get("perfectG") or "").strip() != "1":
            continue
        raw = (row.get("Allele") or "").strip()
        if ";" in raw:
            # Ambiguity list: take the first, but say so.
            first = raw.split(";", 1)[0].strip()
            print(
                f"WARNING: {locus} haplotype has an ambiguous allele list "
                f"({raw}); using {first}",
                file=sys.stderr,
            )
            raw = first
        allele = format_allele(raw)
        if allele is None:
            print(f"WARNING: could not format allele {raw!r}; skipping",
                  file=sys.stderr)
            continue
        by_locus.setdefault(locus, []).append(allele)

    alleles = []
    for locus in CLASS_I_LOCI:
        alleles.extend(by_locus.get(locus, []))

    if not alleles:
        sys.exit(
            f"No confident class I (A/B/C, perfectG=1) alleles in "
            f"{args.bestguess}; cannot run pVACseq"
        )

    print(",".join(alleles))
    print(
        f"Using {len(alleles)} confident class I allele(s): "
        f"{','.join(alleles)}",
        file=sys.stderr,
    )


def cmd_purity(args):
    value = None
    with open_text(args.facets) as handle:
        for line in handle:
            if line.startswith("##purity="):
                value = line.rstrip("\n").split("=", 1)[1].strip()
            if line.startswith("#CHROM"):
                break

    if value is None:
        print("No ##purity= header in FACETS VCF; omitting --tumor-purity",
              file=sys.stderr)
        return
    if value in ("", "NA", ".", "nan", "NaN"):
        print(f"FACETS purity is {value!r}; omitting --tumor-purity",
              file=sys.stderr)
        return

    try:
        purity = float(value)
    except ValueError:
        print(f"FACETS purity {value!r} is not numeric; omitting --tumor-purity",
              file=sys.stderr)
        return

    if not 0 < purity <= 1:
        print(f"FACETS purity {purity} outside (0, 1]; omitting --tumor-purity",
              file=sys.stderr)
        return

    if abs(purity - 0.3) < 1e-10:
        print("FACETS purity is exactly 0.3 (failed-estimate sentinel); "
              "omitting --tumor-purity", file=sys.stderr)
        return

    print(f"--tumor-purity {purity}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare pVACseq inputs for step 8 (VacQ)")
    sub = parser.add_subparsers(dest="command", required=True)

    hla = sub.add_parser("hla", help="confident class I alleles for pVACseq")
    hla.add_argument("bestguess", help="hla/R1_bestguess_G.txt")
    hla.set_defaults(func=cmd_hla)

    purity = sub.add_parser("purity", help="--tumor-purity flag, or empty")
    purity.add_argument("facets", help="FACETS <sample>.facet_out.vcf.gz")
    purity.set_defaults(func=cmd_purity)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
