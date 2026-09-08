#!/usr/bin/env python3
"""Extract KARS G12/Q61 regions and annotate MHC prediction rows.

The program uses only the Python standard library. Region coordinates and the
protein sequence are read directly from SnapGene's Swiss-Prot text export.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import re
import sys
from collections import Counter
from pathlib import Path
from typing import IO, Iterable


SEQUENCE_INPUTS = (
    (1, "V2", "KARS vaccine_V2.sw", "KARS vaccine_V2.fa"),
    (2, "V3", "KARS vaccine_V3.sw", "KARS vaccine_V3.fa"),
    (3, "V4", "KARS vaccine_V4.sw", "KARS vaccine_V4.fa"),
)

REGION_FIELDS = (
    "seq_num",
    "sequence_name",
    "region_name",
    "hotspot",
    "region_start",
    "region_end",
    "mutation_position",
    "mutation_residue",
    "region_sequence",
)

ANNOTATION_FIELDS = (
    "sequence_name",
    "overlaps_G12_Q61_region",
    "within_G12_Q61_region",
    "covers_G12_Q61_mutation",
    "G12_Q61_match_type",
    "G12_Q61_region",
    "G12_Q61_hotspot",
    "G12_Q61_region_start",
    "G12_Q61_region_end",
    "G12_Q61_mutation_position",
    "G12_Q61_covered_mutations",
)

SUMMARY_FIELDS = (
    "seq_num",
    "sequence_name",
    "region_name",
    "hotspot",
    "region_start",
    "region_end",
    "mutation_position",
    "overlap_count",
    "fully_within_count",
    "covers_mutation_count",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract G12/Q61 regions from KARS Swiss-Prot .sw files and add "
            "region markers to an MHC prediction TSV or TSV.GZ."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="Filtered/sorted MHC prediction file (.tsv or .tsv.gz).",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Annotated output path; default is derived from the input name.",
    )
    parser.add_argument(
        "--manifest",
        default="KARS_G12_Q61_regions.tsv",
        help="Extracted region table (default: %(default)s).",
    )
    parser.add_argument(
        "--summary",
        help="Annotation-count summary; default is derived from the output name.",
    )
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Extract and validate the region table without annotating an MHC file.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing derived output files.",
    )
    return parser.parse_args()


def read_single_fasta(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"FASTA file not found: {path}")

    records: list[str] = []
    current: list[str] = []

    with path.open("rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    records.append("".join(current))
                    current = []
                continue
            current.append(line)

    if current:
        records.append("".join(current))

    if len(records) != 1:
        raise ValueError(f"Expected one FASTA record in {path}, found {len(records)}")

    sequence = records[0].rstrip("*").upper()
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYX]+", sequence):
        raise ValueError(f"Invalid amino-acid character in {path}")
    return sequence


def parse_swiss_prot(path: Path) -> tuple[str, list[tuple[str, int, int]]]:
    """Return the protein sequence and labeled G12/Q61 Region features."""
    if not path.is_file():
        raise ValueError(f"Swiss-Prot file not found: {path}")

    sequence_parts: list[str] = []
    labeled_regions: list[tuple[str, int, int]] = []
    current_region: tuple[int, int] | None = None
    in_sequence = False

    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")

            if line.startswith("SQ"):
                in_sequence = True
                current_region = None
                continue

            if in_sequence:
                if line == "//":
                    in_sequence = False
                    continue
                sequence_part = re.sub(r"[\s0-9]", "", line)
                if sequence_part:
                    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYX*]+", sequence_part):
                        raise ValueError(
                            f"Invalid sequence line {line_number} in {path}: {line!r}"
                        )
                    sequence_parts.append(sequence_part)
                continue

            feature_match = re.match(
                r"^FT\s+(\S+)\s+([0-9]+)\s+([0-9]+)\s*$", line
            )
            if feature_match:
                feature_type, start_text, end_text = feature_match.groups()
                current_region = (
                    (int(start_text), int(end_text))
                    if feature_type == "Region"
                    else None
                )
                continue

            label_match = re.match(r'^FT\s+/label="?([^"\s]+)"?\s*$', line)
            if label_match and current_region is not None:
                label = label_match.group(1)
                if re.fullmatch(r"(G12|Q61)[A-Z]", label):
                    labeled_regions.append((label, *current_region))

    if in_sequence:
        raise ValueError(f"Swiss-Prot sequence in {path} has no terminating // line")
    if not sequence_parts:
        raise ValueError(f"No SQ protein sequence found in {path}")
    if not labeled_regions:
        raise ValueError(f"No labeled G12/Q61 Region features found in {path}")

    return "".join(sequence_parts).rstrip("*"), labeled_regions


def extract_regions(
    seq_num: int,
    sequence_name: str,
    swiss_prot_path: Path,
    fasta_path: Path,
) -> list[dict[str, object]]:
    fasta_sequence = read_single_fasta(fasta_path)
    swiss_prot_sequence, labeled_regions = parse_swiss_prot(swiss_prot_path)

    if swiss_prot_sequence != fasta_sequence:
        raise ValueError(
            f"Protein sequence mismatch between {swiss_prot_path} and {fasta_path}"
        )

    regions: list[dict[str, object]] = []

    for region_name, region_start, region_end in labeled_regions:
        match = re.fullmatch(r"(G12|Q61)([A-Z])", region_name)
        if match is None:
            raise ValueError(
                f"Unexpected region label {region_name!r} in {swiss_prot_path}"
            )
        hotspot, expected_residue = match.groups()
        mutation_offset = 11 if hotspot == "G12" else 6
        expected_length = 25 if hotspot == "G12" else 27
        mutation_position = region_start + mutation_offset

        if region_end - region_start + 1 != expected_length:
            raise ValueError(
                f"Unexpected {region_name} region length at "
                f"{region_start}-{region_end} in {swiss_prot_path}"
            )
        if region_end > len(fasta_sequence):
            raise ValueError(
                f"Region {region_start}-{region_end} exceeds {fasta_path} length"
            )

        actual_residue = fasta_sequence[mutation_position - 1]
        if actual_residue != expected_residue:
            raise ValueError(
                f"Residue validation failed for {sequence_name} {region_name}: "
                f"position {mutation_position} is {actual_residue}, expected "
                f"{expected_residue}"
            )

        regions.append(
            {
                "seq_num": seq_num,
                "sequence_name": sequence_name,
                "region_name": region_name,
                "hotspot": hotspot,
                "region_start": region_start,
                "region_end": region_end,
                "mutation_position": mutation_position,
                "mutation_residue": actual_residue,
                "region_sequence": fasta_sequence[region_start - 1 : region_end],
            }
        )

    regions.sort(key=lambda row: (int(row["region_start"]), int(row["region_end"])))
    if not regions:
        raise ValueError(f"No G12/Q61 regions found in {swiss_prot_path}")

    for previous, current in zip(regions, regions[1:]):
        if int(current["region_start"]) <= int(previous["region_end"]):
            raise ValueError(
                f"Overlapping regions in {swiss_prot_path}: "
                f"{previous['region_name']} and {current['region_name']}"
            )

    return regions


def build_region_table(base_dir: Path) -> list[dict[str, object]]:
    regions: list[dict[str, object]] = []
    for seq_num, sequence_name, swiss_prot_name, fasta_name in SEQUENCE_INPUTS:
        regions.extend(
            extract_regions(
                seq_num,
                sequence_name,
                base_dir / swiss_prot_name,
                base_dir / fasta_name,
            )
        )
    return regions


def temporary_path(destination: Path) -> Path:
    return destination.with_name(f".{destination.name}.partial.{os.getpid()}")


def ensure_writable_destination(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise ValueError(f"Output already exists (use --force to replace it): {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def write_tsv_atomic(
    path: Path,
    fields: Iterable[str],
    rows: Iterable[dict[str, object]],
    force: bool,
    allow_identical: bool = False,
) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(fields),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    content = buffer.getvalue()

    if path.exists() and not force:
        if allow_identical and path.read_text(encoding="utf-8") == content:
            return
        raise ValueError(f"Output already exists (use --force to replace it): {path}")

    path.parent.mkdir(parents=True, exist_ok=True)
    partial = temporary_path(path)
    try:
        with partial.open("wt", encoding="utf-8", newline="") as handle:
            handle.write(content)
        os.replace(partial, path)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise


def derive_annotated_output(input_path: Path) -> Path:
    name = input_path.name
    if name.endswith(".tsv.gz"):
        return input_path.with_name(
            f"{name[:-7]}_G12_Q61_annotated.tsv.gz"
        )
    if name.endswith(".tsv"):
        return input_path.with_name(f"{name[:-4]}_G12_Q61_annotated.tsv")
    if name.endswith(".gz"):
        return input_path.with_name(f"{name[:-3]}_G12_Q61_annotated.gz")
    return input_path.with_name(f"{name}_G12_Q61_annotated.tsv")


def derive_summary_output(output_path: Path) -> Path:
    name = output_path.name
    if name.endswith(".tsv.gz"):
        name = name[:-7]
    elif name.endswith(".tsv"):
        name = name[:-4]
    elif name.endswith(".gz"):
        name = name[:-3]
    return output_path.with_name(f"{name}_summary.tsv")


def open_input_text(path: Path) -> IO[str]:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def open_output_text(path: Path, compressed: bool) -> IO[str]:
    if compressed:
        return gzip.open(path, "wt", encoding="utf-8", newline="", compresslevel=1)
    return path.open("wt", encoding="utf-8", newline="")


def integer_field(row: dict[str, str], field: str, row_number: int) -> int:
    try:
        return int(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid {field!r} value on input row {row_number}: {row.get(field)!r}"
        ) from exc


def annotate_predictions(
    input_path: Path,
    output_path: Path,
    summary_path: Path,
    regions: list[dict[str, object]],
    force: bool,
) -> tuple[int, int, int, int]:
    if not input_path.is_file():
        raise ValueError(f"MHC input file not found: {input_path}")

    ensure_writable_destination(output_path, force)
    ensure_writable_destination(summary_path, force)

    regions_by_seq: dict[int, list[dict[str, object]]] = {}
    sequence_names: dict[int, str] = {}
    for region in regions:
        seq_num = int(region["seq_num"])
        regions_by_seq.setdefault(seq_num, []).append(region)
        sequence_names[seq_num] = str(region["sequence_name"])

    counts_within: Counter[tuple[int, int, int, str]] = Counter()
    counts_overlap: Counter[tuple[int, int, int, str]] = Counter()
    counts_covering: Counter[tuple[int, int, int, str]] = Counter()
    total_rows = 0
    total_overlap = 0
    total_within = 0
    total_covering = 0

    output_partial = temporary_path(output_path)
    summary_partial = temporary_path(summary_path)
    compressed_output = output_path.name.endswith(".gz")

    try:
        with open_input_text(input_path) as source, open_output_text(
            output_partial, compressed_output
        ) as destination:
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"MHC input has no header: {input_path}")

            required = {"seq_num", "start", "end"}
            missing = required.difference(reader.fieldnames)
            if missing:
                raise ValueError(
                    f"MHC input is missing required columns: {', '.join(sorted(missing))}"
                )

            already_present = set(ANNOTATION_FIELDS).intersection(reader.fieldnames)
            if already_present:
                raise ValueError(
                    "MHC input already contains annotation columns: "
                    + ", ".join(sorted(already_present))
                )

            output_fields = list(reader.fieldnames) + list(ANNOTATION_FIELDS)
            writer = csv.DictWriter(
                destination,
                fieldnames=output_fields,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()

            for row_number, row in enumerate(reader, start=2):
                total_rows += 1
                seq_num = integer_field(row, "seq_num", row_number)
                peptide_start = integer_field(row, "start", row_number)
                peptide_end = integer_field(row, "end", row_number)

                if seq_num not in regions_by_seq:
                    raise ValueError(
                        f"Unexpected seq_num {seq_num} on input row {row_number}; "
                        "expected 1=V2, 2=V3, or 3=V4"
                    )
                if peptide_start < 1 or peptide_end < peptide_start:
                    raise ValueError(
                        f"Invalid peptide interval {peptide_start}-{peptide_end} "
                        f"on input row {row_number}"
                    )

                if "length" in row and row["length"] not in (None, ""):
                    peptide_length = integer_field(row, "length", row_number)
                    if peptide_end - peptide_start + 1 != peptide_length:
                        raise ValueError(
                            f"Coordinate/length mismatch on input row {row_number}: "
                            f"{peptide_start}-{peptide_end}, length={peptide_length}"
                        )

                overlapping_regions = [
                    region
                    for region in regions_by_seq[seq_num]
                    if peptide_start <= int(region["region_end"])
                    and peptide_end >= int(region["region_start"])
                ]

                fully_within_regions = [
                    region
                    for region in overlapping_regions
                    if peptide_start >= int(region["region_start"])
                    and peptide_end <= int(region["region_end"])
                ]

                covered_mutation_regions = [
                    region
                    for region in overlapping_regions
                    if peptide_start
                    <= int(region["mutation_position"])
                    <= peptide_end
                ]

                has_overlap = bool(overlapping_regions)
                is_fully_within = bool(fully_within_regions)
                covers_mutation = bool(covered_mutation_regions)

                if is_fully_within and covers_mutation:
                    match_type = "fully_within_and_covers"
                elif is_fully_within:
                    match_type = "fully_within_not_cover"
                elif has_overlap and covers_mutation:
                    match_type = "partial_overlap_covers"
                elif has_overlap:
                    match_type = "partial_overlap_not_cover"
                else:
                    match_type = "no_overlap"

                row["sequence_name"] = sequence_names[seq_num]
                row["overlaps_G12_Q61_region"] = "1" if has_overlap else "0"
                row["within_G12_Q61_region"] = "1" if is_fully_within else "0"
                row["covers_G12_Q61_mutation"] = "1" if covers_mutation else "0"
                row["G12_Q61_match_type"] = match_type

                if has_overlap:
                    row["G12_Q61_region"] = ";".join(
                        str(region["region_name"]) for region in overlapping_regions
                    )
                    row["G12_Q61_hotspot"] = ";".join(
                        str(region["hotspot"]) for region in overlapping_regions
                    )
                    row["G12_Q61_region_start"] = ";".join(
                        str(region["region_start"]) for region in overlapping_regions
                    )
                    row["G12_Q61_region_end"] = ";".join(
                        str(region["region_end"]) for region in overlapping_regions
                    )
                    row["G12_Q61_mutation_position"] = ";".join(
                        str(region["mutation_position"])
                        for region in overlapping_regions
                    )
                    row["G12_Q61_covered_mutations"] = ";".join(
                        f"{region['region_name']}@{region['mutation_position']}"
                        for region in covered_mutation_regions
                    )

                    total_overlap += 1

                for region in overlapping_regions:
                    region_key = (
                        seq_num,
                        int(region["region_start"]),
                        int(region["region_end"]),
                        str(region["region_name"]),
                    )
                    counts_overlap[region_key] += 1

                for region in fully_within_regions:
                    region_key = (
                        seq_num,
                        int(region["region_start"]),
                        int(region["region_end"]),
                        str(region["region_name"]),
                    )
                    counts_within[region_key] += 1
                if is_fully_within:
                    total_within += 1

                for region in covered_mutation_regions:
                    region_key = (
                        seq_num,
                        int(region["region_start"]),
                        int(region["region_end"]),
                        str(region["region_name"]),
                    )
                    counts_covering[region_key] += 1
                if covers_mutation:
                    total_covering += 1

                if not has_overlap:
                    row["G12_Q61_region"] = ""
                    row["G12_Q61_hotspot"] = ""
                    row["G12_Q61_region_start"] = ""
                    row["G12_Q61_region_end"] = ""
                    row["G12_Q61_mutation_position"] = ""
                    row["G12_Q61_covered_mutations"] = ""

                writer.writerow(row)

        summary_rows: list[dict[str, object]] = [
            {
                "seq_num": "ALL",
                "sequence_name": "ALL",
                "region_name": "ALL",
                "hotspot": "ALL",
                "region_start": "",
                "region_end": "",
                "mutation_position": "",
                "overlap_count": total_overlap,
                "fully_within_count": total_within,
                "covers_mutation_count": total_covering,
            }
        ]

        for region in regions:
            region_key = (
                int(region["seq_num"]),
                int(region["region_start"]),
                int(region["region_end"]),
                str(region["region_name"]),
            )
            summary_rows.append(
                {
                    "seq_num": region["seq_num"],
                    "sequence_name": region["sequence_name"],
                    "region_name": region["region_name"],
                    "hotspot": region["hotspot"],
                    "region_start": region["region_start"],
                    "region_end": region["region_end"],
                    "mutation_position": region["mutation_position"],
                    "overlap_count": counts_overlap[region_key],
                    "fully_within_count": counts_within[region_key],
                    "covers_mutation_count": counts_covering[region_key],
                }
            )

        with summary_partial.open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(SUMMARY_FIELDS),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(summary_rows)

        os.replace(output_partial, output_path)
        os.replace(summary_partial, summary_path)
    except BaseException:
        output_partial.unlink(missing_ok=True)
        summary_partial.unlink(missing_ok=True)
        raise

    return total_rows, total_overlap, total_within, total_covering


def main() -> int:
    args = parse_args()
    base_dir = Path.cwd()

    try:
        regions = build_region_table(base_dir)
        manifest_path = Path(args.manifest)
        write_tsv_atomic(
            manifest_path,
            REGION_FIELDS,
            regions,
            force=args.force,
            allow_identical=True,
        )

        counts_by_sequence = Counter(str(row["sequence_name"]) for row in regions)
        print(
            "Extracted regions: "
            + ", ".join(
                f"{name}={counts_by_sequence[name]}" for name in ("V2", "V3", "V4")
            ),
            file=sys.stderr,
        )
        print(f"Region manifest: {manifest_path.resolve()}", file=sys.stderr)

        if args.manifest_only:
            return 0
        if not args.input:
            raise ValueError("MHC input is required unless --manifest-only is used")

        input_path = Path(args.input)
        output_path = Path(args.output) if args.output else derive_annotated_output(input_path)
        summary_path = (
            Path(args.summary) if args.summary else derive_summary_output(output_path)
        )

        total, overlapping, within, covering = annotate_predictions(
            input_path,
            output_path,
            summary_path,
            regions,
            force=args.force,
        )
        print(f"Rows processed: {total}", file=sys.stderr)
        print(f"Rows overlapping a G12/Q61 region: {overlapping}", file=sys.stderr)
        print(f"Rows fully within a G12/Q61 region: {within}", file=sys.stderr)
        print(f"Rows also covering the mutation position: {covering}", file=sys.stderr)
        print(f"Annotated result: {output_path.resolve()}", file=sys.stderr)
        print(f"Annotation summary: {summary_path.resolve()}", file=sys.stderr)
        return 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
