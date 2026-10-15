#!/usr/bin/env python3
"""Score and rank pVACseq MHC-I 9-mer neoantigen candidates.

The script implements two published, distinct quality definitions:

1. NeoantigenEditing 2022 quality:
   Q22 = R * [w * dC(MT, WT) + (1-w) * ln(KdWT/KdMT)]

2. 2023 vaccine-response 9-mer quality:
   Q23(MT9) = dC(MT9, WT9) - min_{IEDB9} dC(MT9, IEDB9)

It reads pVACseq ``MHC_I.all_epitopes.tsv`` output, collapses duplicate
transcript rows, optionally merges RNA and clonality/CCF tables, and writes
machine-readable TSV reports plus a NeoantigenEditing peptide JSON.

Only Python's standard library and local NCBI BLAST+ executables are required.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


AA = set("ACDEFGHIKLMNPQRSTVWY")
MISSING = {"", "NA", "N/A", "NONE", "NULL", "."}

# BLOSUM62, used by the original NeoantigenEditing alignment code.
_BLOSUM62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -3 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""


def parse_blosum62() -> dict[tuple[str, str], int]:
    lines = [line.split() for line in _BLOSUM62_TEXT.strip().splitlines()]
    header = lines[0]
    result: dict[tuple[str, str], int] = {}
    for row in lines[1:]:
        aa = row[0]
        for bb, value in zip(header, row[1:]):
            result[(aa, bb)] = int(value)
    return result


BLOSUM62 = parse_blosum62()


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    return str(value).strip().upper() in MISSING


def to_float(value: Any) -> float | None:
    if is_missing(value):
        return None
    try:
        x = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def to_int(value: Any) -> int | None:
    x = to_float(value)
    return int(x) if x is not None else None


def clean_text(value: Any) -> str:
    return "" if is_missing(value) else str(value).strip()


def norm_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "clonal"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            out: dict[str, Any] = {}
            for field in fields:
                value = row.get(field)
                if value is None:
                    out[field] = "NA"
                elif isinstance(value, float):
                    out[field] = f"{value:.8g}"
                else:
                    out[field] = value
            writer.writerow(out)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(chunks).upper()
                identifier = line[1:].split()[0].split("|")[0]
                chunks = []
            else:
                chunks.append(line)
    if identifier is not None:
        records[identifier] = "".join(chunks).upper()
    return records


def valid_peptide(seq: str, length: int = 9) -> bool:
    return len(seq) == length and set(seq) <= AA


class EpitopeDistance:
    def __init__(self, model_json: Path):
        with model_json.open("r", encoding="utf-8") as handle:
            model = json.load(handle)
        self.d_i = [float(x) for x in model["d_i"]]
        self.matrix = {key: float(value) for key, value in model["M_ab"].items()}
        if len(self.d_i) != 9:
            raise ValueError("Epitope-distance model must contain 9 position weights")

    def distance(self, seq_a: str, seq_b: str) -> float:
        if not valid_peptide(seq_a) or not valid_peptide(seq_b):
            raise ValueError("Epitope distance requires two canonical 9-mer peptides")
        return sum(
            self.d_i[i] * self.matrix[f"{seq_a[i]}->{seq_b[i]}"]
            for i in range(9)
        )


def local_alignment_score(
    seq_a: str,
    seq_b: str,
    gap_open: int = -11,
    gap_extend: int = -1,
) -> float:
    """Smith-Waterman local alignment score with affine gaps."""
    n, m = len(seq_a), len(seq_b)
    h_prev = [0.0] * (m + 1)
    gap_a_prev = [float("-inf")] * (m + 1)
    best = 0.0
    for i in range(1, n + 1):
        h_cur = [0.0] * (m + 1)
        gap_a_cur = [float("-inf")] * (m + 1)
        gap_b = float("-inf")
        for j in range(1, m + 1):
            gap_a_cur[j] = max(
                h_prev[j] + gap_open,
                gap_a_prev[j] + gap_extend,
            )
            gap_b = max(h_cur[j - 1] + gap_open, gap_b + gap_extend)
            subst = BLOSUM62[(seq_a[i - 1], seq_b[j - 1])]
            h_cur[j] = max(0.0, h_prev[j - 1] + subst, gap_a_cur[j], gap_b)
            best = max(best, h_cur[j])
        h_prev = h_cur
        gap_a_prev = gap_a_cur
    return best


def blast_candidates(
    peptides: Iterable[str],
    iedb_fasta: Path,
    iedb_records: dict[str, str],
    output_dir: Path,
    blastp_exe: str,
    makeblastdb_exe: str,
) -> dict[str, list[float]]:
    peptide_list = sorted(set(peptides))
    if not peptide_list:
        return {}
    if shutil.which(blastp_exe) is None and not Path(blastp_exe).exists():
        raise FileNotFoundError(f"blastp not found: {blastp_exe}")

    existing_db = all(Path(str(iedb_fasta) + suffix).exists() for suffix in (".pin", ".phr", ".psq"))

    with tempfile.TemporaryDirectory(prefix="neoedit_", dir=output_dir) as temp_name:
        temp_dir = Path(temp_name)
        query_path = temp_dir / "query.fasta"
        result_path = temp_dir / "blast.tsv"
        query_to_peptide: dict[str, str] = {}
        with query_path.open("w", encoding="utf-8") as handle:
            for i, peptide in enumerate(peptide_list):
                query_id = f"q{i}"
                query_to_peptide[query_id] = peptide
                handle.write(f">{query_id}\n{peptide}\n")

        db_path = iedb_fasta
        if not existing_db:
            if shutil.which(makeblastdb_exe) is None and not Path(makeblastdb_exe).exists():
                raise FileNotFoundError(f"makeblastdb not found: {makeblastdb_exe}")
            db_path = temp_dir / "iedb_db"
            subprocess.run(
                [
                    makeblastdb_exe,
                    "-in",
                    str(iedb_fasta),
                    "-dbtype",
                    "prot",
                    "-out",
                    str(db_path),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )

        cmd = [
            blastp_exe,
            "-db",
            str(db_path),
            "-query",
            str(query_path),
            "-outfmt",
            "6 qseqid sacc score",
            "-gapopen",
            "32767",
            "-gapextend",
            "32767",
            "-evalue",
            "1e6",
            "-max_hsps",
            "1",
            "-matrix",
            "BLOSUM62",
            "-max_target_seqs",
            "10000000",
            "-out",
            str(result_path),
        ]
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)

        peptide_to_ids: dict[str, set[str]] = defaultdict(set)
        with result_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2 or parts[0] not in query_to_peptide:
                    continue
                epitope_id = parts[1].split("|")[0]
                peptide_to_ids[query_to_peptide[parts[0]]].add(epitope_id)

    scores: dict[str, list[float]] = {}
    for peptide in peptide_list:
        values: list[float] = []
        for epitope_id in peptide_to_ids.get(peptide, set()):
            epitope = iedb_records.get(epitope_id)
            if not epitope or not set(epitope) <= AA:
                continue
            values.append(local_alignment_score(peptide, epitope))
        scores[peptide] = values
    return scores


def recognition_potential(scores: list[float], a: float = 22.897590714815188, k: float = 1.0) -> float:
    if not scores:
        return 0.0
    values = [-k * (a - score) for score in scores]
    maximum = max(values)
    log_b = maximum + math.log(sum(math.exp(x - maximum) for x in values))
    if log_b >= 0:
        return 1.0 / (1.0 + math.exp(-log_b))
    exp_log_b = math.exp(log_b)
    return exp_log_b / (1.0 + exp_log_b)


def nearest_iedb_9mer(
    peptide: str,
    iedb_9mers: list[str],
    distance_model: EpitopeDistance,
) -> tuple[float, str]:
    best_distance = float("inf")
    best_sequence = ""
    for target in iedb_9mers:
        value = distance_model.distance(peptide, target)
        if value < best_distance:
            best_distance = value
            best_sequence = target
    return best_distance, best_sequence


def first_existing(row: dict[str, Any], names: Iterable[str]) -> Any:
    normalized = {norm_header(k): v for k, v in row.items()}
    for name in names:
        value = normalized.get(norm_header(name))
        if not is_missing(value):
            return value
    return None


def transcript_priority(row: dict[str, str]) -> tuple[int, int, int, int]:
    mane = 1 if clean_text(row.get("MANE Select")) not in {"", "False"} else 0
    canonical = 1 if truthy(row.get("Canonical")) else 0
    tsl = to_int(row.get("Transcript Support Level"))
    transcript_length = to_int(row.get("Transcript Length")) or 0
    return mane, canonical, -(tsl if tsl is not None else 999), transcript_length


def pvac_variant_id(row: dict[str, Any]) -> str:
    return "-".join(
        clean_text(row.get(field))
        for field in ("Chromosome", "Start", "Stop", "Reference", "Variant")
    )


def mutation_id(row: dict[str, Any]) -> str:
    chrom = clean_text(row.get("Chromosome"))
    start = to_int(row.get("Start"))
    position = start + 1 if start is not None else clean_text(row.get("Start"))
    return f"{chrom}:{position}:{clean_text(row.get('Reference'))}>{clean_text(row.get('Variant'))}"


def load_aggregated(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    result: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        key = clean_text(row.get("ID"))
        if key:
            result[key] = row
    return result


def load_filtered_candidate_keys(path: Path | None) -> set[tuple[str, str, str]]:
    if path is None or not path.exists():
        return set()
    keys: set[tuple[str, str, str]] = set()
    for row in read_tsv(path):
        keys.add(
            (
                pvac_variant_id(row),
                clean_text(row.get("HLA Allele")),
                clean_text(row.get("MT Epitope Seq")).upper(),
            )
        )
    return keys


def load_supplement(path: Path | None) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    by_mutation: dict[str, dict[str, str]] = {}
    by_gene: dict[str, dict[str, str]] = {}
    if path is None:
        return by_mutation, by_gene
    for raw in read_tsv(path):
        row = {norm_header(k): v for k, v in raw.items()}
        key = first_existing(row, ["mutation_id", "pvac_variant_id", "variant_id"])
        gene = first_existing(row, ["gene", "gene_name", "symbol"])
        if key:
            by_mutation[str(key)] = row
        if gene:
            by_gene[str(gene)] = row
    return by_mutation, by_gene


def supplement_for(
    candidate: dict[str, Any],
    by_mutation: dict[str, dict[str, str]],
    by_gene: dict[str, dict[str, str]],
) -> dict[str, str]:
    for key in (candidate["mutation_id"], candidate["pvac_variant_id"]):
        if key in by_mutation:
            return by_mutation[key]
    return by_gene.get(candidate["gene"], {})


def percentile_scores(values: list[float | None]) -> list[float | None]:
    indexed = sorted((value, i) for i, value in enumerate(values) if value is not None)
    output: list[float | None] = [None] * len(values)
    if not indexed:
        return output
    if len(indexed) == 1:
        output[indexed[0][1]] = 1.0
        return output
    pos = 0
    while pos < len(indexed):
        end = pos + 1
        while end < len(indexed) and indexed[end][0] == indexed[pos][0]:
            end += 1
        average_rank = ((pos + 1) + end) / 2.0
        score = (average_rank - 1.0) / (len(indexed) - 1.0)
        for _, original_index in indexed[pos:end]:
            output[original_index] = score
        pos = end
    return output


def parse_weights(spec: str) -> dict[str, float]:
    expected = {"q2023", "q2022", "binding", "rna_expr", "rna_vaf", "dna_vaf", "clonality"}
    result: dict[str, float] = {}
    for item in spec.split(","):
        key, sep, value = item.partition("=")
        if not sep or key.strip() not in expected:
            raise ValueError(f"Invalid weight entry: {item}")
        result[key.strip()] = float(value)
    missing = expected - result.keys()
    if missing:
        raise ValueError(f"Missing weights: {', '.join(sorted(missing))}")
    if any(value < 0 for value in result.values()) or sum(result.values()) <= 0:
        raise ValueError("Weights must be non-negative and have a positive sum")
    total = sum(result.values())
    return {key: value / total for key, value in result.items()}


def build_candidates(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    required = {
        "Chromosome",
        "Start",
        "Stop",
        "Reference",
        "Variant",
        "HLA Allele",
        "MT Epitope Seq",
        "WT Epitope Seq",
        "Best MT IC50 Score",
        "Corresponding WT IC50 Score",
    }
    if not rows:
        raise ValueError("pVACseq input has no data rows")
    absent = sorted(required - rows[0].keys())
    if absent:
        raise ValueError(f"Missing required pVACseq columns: {', '.join(absent)}")

    selected: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        mt = clean_text(row.get("MT Epitope Seq")).upper()
        wt = clean_text(row.get("WT Epitope Seq")).upper()
        # Frameshift peptides usually have no meaningful position-matched WT
        # peptide in pVACseq.  Keep them for binding/ranking, but leave the two
        # published MT-vs-WT quality scores unavailable.
        if not valid_peptide(mt):
            continue
        key = (
            clean_text(row.get("Chromosome")),
            clean_text(row.get("Start")),
            clean_text(row.get("Stop")),
            clean_text(row.get("Reference")),
            clean_text(row.get("Variant")),
            clean_text(row.get("HLA Allele")),
            mt,
            wt if valid_peptide(wt) else "",
        )
        if key not in selected or transcript_priority(row) > transcript_priority(selected[key]):
            selected[key] = row

    candidates: list[dict[str, Any]] = []
    for row in selected.values():
        mt = clean_text(row.get("MT Epitope Seq")).upper()
        wt = clean_text(row.get("WT Epitope Seq")).upper()
        if not valid_peptide(wt):
            wt = ""
        mt_ic50 = to_float(row.get("NetMHCpan MT IC50 Score"))
        if mt_ic50 is None:
            mt_ic50 = to_float(row.get("Best MT IC50 Score"))
        wt_ic50 = to_float(row.get("NetMHCpan WT IC50 Score"))
        if wt_ic50 is None:
            wt_ic50 = to_float(row.get("Corresponding WT IC50 Score"))
        mt_percentile = to_float(row.get("NetMHCpan MT Percentile"))
        if mt_percentile is None:
            mt_percentile = to_float(row.get("Best MT Percentile"))
        start = to_int(row.get("Start"))
        candidate = {
            "sample_id": "",
            "mutation_id": mutation_id(row),
            "pvac_variant_id": pvac_variant_id(row),
            "chromosome": clean_text(row.get("Chromosome")),
            "start_0based": start,
            "position_1based": start + 1 if start is not None else None,
            "stop_0based": to_int(row.get("Stop")),
            "reference": clean_text(row.get("Reference")),
            "variant": clean_text(row.get("Variant")),
            "variant_type": clean_text(row.get("Variant Type")),
            "gene": clean_text(row.get("Gene Name")),
            "ensembl_gene_id": clean_text(row.get("Ensembl Gene ID")),
            "transcript": clean_text(row.get("Transcript")),
            "transcript_support_level": to_int(row.get("Transcript Support Level")),
            "canonical": clean_text(row.get("Canonical")),
            "mane_select": clean_text(row.get("MANE Select")),
            "hgvsc": clean_text(row.get("HGVSc")),
            "hgvsp": clean_text(row.get("HGVSp")),
            "amino_acid_change": clean_text(row.get("Mutation")),
            "protein_position": clean_text(row.get("Protein Position")),
            "hla_allele": clean_text(row.get("HLA Allele")),
            "mt_peptide": mt,
            "wt_peptide": wt,
            "mutation_position_in_peptide": to_int(row.get("Mutation Position")),
            "mt_ic50_nm": mt_ic50,
            "wt_ic50_nm": wt_ic50,
            "mt_percentile": mt_percentile,
            "wt_percentile": to_float(row.get("NetMHCpan WT Percentile"))
            or to_float(row.get("Corresponding WT Percentile")),
            "tumor_dna_depth": to_int(row.get("Tumor DNA Depth")),
            "tumor_dna_vaf": to_float(row.get("Tumor DNA VAF")),
            "tumor_rna_depth": to_int(row.get("Tumor RNA Depth")),
            "tumor_rna_vaf": to_float(row.get("Tumor RNA VAF")),
            "normal_depth": to_int(row.get("Normal Depth")),
            "normal_vaf": to_float(row.get("Normal VAF")),
            "gene_expression": to_float(row.get("Gene Expression")),
            "transcript_expression": to_float(row.get("Transcript Expression")),
        }
        candidates.append(candidate)
    return candidates


def apply_supplements(
    candidates: list[dict[str, Any]],
    rna_maps: tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]],
    clone_maps: tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]],
    clonal_ccf_threshold: float,
) -> None:
    rna_by_mutation, rna_by_gene = rna_maps
    clone_by_mutation, clone_by_gene = clone_maps
    for candidate in candidates:
        rna = supplement_for(candidate, rna_by_mutation, rna_by_gene)
        clone = supplement_for(candidate, clone_by_mutation, clone_by_gene)

        for output_name, aliases in {
            "gene_expression": ["gene_expression", "gene_expr", "tpm", "fpkm"],
            "transcript_expression": ["transcript_expression", "transcript_tpm", "transcript_fpkm"],
            "tumor_rna_depth": ["rna_depth", "tumor_rna_depth", "trna_depth"],
            "tumor_rna_vaf": ["rna_vaf", "tumor_rna_vaf", "trna_vaf"],
        }.items():
            value = to_float(first_existing(rna, aliases))
            if value is not None:
                candidate[output_name] = int(value) if output_name.endswith("depth") else value

        ccf = to_float(first_existing(clone, ["ccf", "cancer_cell_fraction", "cellular_prevalence"]))
        clone_id = first_existing(clone, ["clone_id", "cluster_id", "cluster"])
        clonality = first_existing(clone, ["clonality", "clonal_status", "is_clonal"])
        candidate["ccf"] = ccf
        candidate["clone_id"] = clean_text(clone_id)
        candidate["best_tree_cellular_prevalence"] = to_float(
            first_existing(clone, ["best_tree_cellular_prevalence", "cellular_prevalence"])
        )
        candidate["best_tree_ccf"] = to_float(first_existing(clone, ["best_tree_ccf"]))
        candidate["ccf_min_top_trees"] = to_float(first_existing(clone, ["ccf_min_top_trees"]))
        candidate["ccf_max_top_trees"] = to_float(first_existing(clone, ["ccf_max_top_trees"]))
        candidate["n_trees_matched"] = to_int(first_existing(clone, ["n_trees_matched"]))
        candidate["tree_assignments"] = clean_text(first_existing(clone, ["tree_assignments"]))
        if clonality is not None:
            if truthy(clonality):
                candidate["clonality"] = "clonal"
            elif str(clonality).strip().lower() in {"0", "false", "no", "n", "subclonal"}:
                candidate["clonality"] = "subclonal"
            else:
                candidate["clonality"] = clean_text(clonality)
            candidate["clonality_source"] = "input"
        elif ccf is not None:
            candidate["clonality"] = "clonal" if ccf >= clonal_ccf_threshold else "subclonal"
            candidate["clonality_source"] = f"derived_ccf_threshold_{clonal_ccf_threshold:g}"
        else:
            candidate["clonality"] = ""
            candidate["clonality_source"] = ""


def add_quality_scores(
    candidates: list[dict[str, Any]],
    distance_model: EpitopeDistance,
    iedb_9mers: list[str],
    blast_scores: dict[str, list[float]],
    w: float = 0.22402192838740312,
) -> None:
    q23_cache: dict[str, tuple[float, str]] = {}
    for candidate in candidates:
        mt = candidate["mt_peptide"]
        wt = candidate["wt_peptide"]
        dc_mt_wt = distance_model.distance(mt, wt) if valid_peptide(wt) else None
        if mt not in q23_cache:
            q23_cache[mt] = nearest_iedb_9mer(mt, iedb_9mers, distance_model)
        min_iedb_distance, nearest_iedb = q23_cache[mt]
        r_value = recognition_potential(blast_scores.get(mt, []))
        mt_ic50 = candidate["mt_ic50_nm"]
        wt_ic50 = candidate["wt_ic50_nm"]
        log_a = None
        q22 = None
        if (
            dc_mt_wt is not None
            and mt_ic50 is not None
            and wt_ic50 is not None
            and mt_ic50 > 0
            and wt_ic50 > 0
        ):
            log_a = math.log(wt_ic50 / mt_ic50)
            q22 = (w * dc_mt_wt + (1.0 - w) * log_a) * r_value
        candidate.update(
            {
                "q2022_R": r_value,
                "q2022_logC": dc_mt_wt,
                "q2022_logA": log_a,
                "quality_2022": q22,
                "q2023_dc_mt_wt": dc_mt_wt,
                "q2023_min_dc_iedb": min_iedb_distance,
                "q2023_nearest_iedb_9mer": nearest_iedb,
                "quality_2023_9mer": (
                    dc_mt_wt - min_iedb_distance if dc_mt_wt is not None else None
                ),
            }
        )


def add_eligibility_and_flags(
    candidates: list[dict[str, Any]],
    q22_ic50_cutoff: float,
    q23_ic50_cutoff: float,
    pvac_binding_ic50_cutoff: float,
    binding_percentile_cutoff: float,
    min_tdna_depth: int,
    min_tdna_vaf: float,
    min_normal_depth: int,
    max_normal_vaf: float,
) -> None:
    for row in candidates:
        ic50 = row["mt_ic50_nm"]
        percentile = row["mt_percentile"]
        row["eligible_2022_ic50"] = bool(ic50 is not None and ic50 <= q22_ic50_cutoff)
        row["eligible_2023_ic50"] = bool(ic50 is not None and ic50 <= q23_ic50_cutoff)
        row["pvac_binding_pass"] = bool(
            ic50 is not None
            and percentile is not None
            and ic50 <= pvac_binding_ic50_cutoff
            and percentile <= binding_percentile_cutoff
        )
        row["analysis_tdna_support_pass"] = bool(
            row["tumor_dna_depth"] is not None
            and row["tumor_dna_vaf"] is not None
            and row["tumor_dna_depth"] >= min_tdna_depth
            and row["tumor_dna_vaf"] >= min_tdna_vaf
        )
        row["normal_support_pass"] = bool(
            row["normal_depth"] is not None
            and row["normal_vaf"] is not None
            and row["normal_depth"] >= min_normal_depth
            and row["normal_vaf"] < max_normal_vaf
        )
        flags: list[str] = []
        if not row["analysis_tdna_support_pass"]:
            flags.append("low_tumor_DNA_support")
        if not row["normal_support_pass"]:
            flags.append("normal_evidence_review")
        if row["gene_expression"] is None and row["transcript_expression"] is None:
            flags.append("RNA_expression_missing")
        if row["tumor_rna_vaf"] is None:
            flags.append("RNA_VAF_missing")
        if row.get("ccf") is None and not row.get("clonality"):
            flags.append("CCF_clonality_missing")
        if (
            row.get("ccf_min_top_trees") is not None
            and row.get("ccf_max_top_trees") is not None
            and row["ccf_min_top_trees"] < 0.9 <= row["ccf_max_top_trees"]
        ):
            flags.append("clonality_top_trees_uncertain")
        if not valid_peptide(row["wt_peptide"]):
            flags.append("quality_unavailable_missing_WT")
        row["review_flags"] = ";".join(flags) if flags else "none"


def add_mutation_level_scores(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for candidate in candidates:
        groups[candidate["mutation_id"]].append(candidate)

    summaries: list[dict[str, Any]] = []
    for key, group in groups.items():
        eligible_23 = [
            row
            for row in group
            if row["eligible_2023_ic50"] and row["quality_2023_9mer"] is not None
        ]
        best_per_peptide: dict[str, dict[str, Any]] = {}
        for row in eligible_23:
            peptide = row["mt_peptide"]
            old = best_per_peptide.get(peptide)
            if old is None or row["mt_ic50_nm"] < old["mt_ic50_nm"]:
                best_per_peptide[peptide] = row
        top_two = sorted(
            best_per_peptide.values(),
            key=lambda row: (row["quality_2023_9mer"], -row["mt_ic50_nm"]),
            reverse=True,
        )[:2]
        q23_mean = (
            sum(row["quality_2023_9mer"] for row in top_two) / len(top_two)
            if top_two
            else None
        )
        q22_eligible = [
            row["quality_2022"]
            for row in group
            if row["eligible_2022_ic50"] and row["quality_2022"] is not None
        ]
        best_binding = min(
            (row for row in group if row["mt_ic50_nm"] is not None),
            key=lambda row: row["mt_ic50_nm"],
        )
        summary = {
            "sample_id": group[0]["sample_id"],
            "mutation_id": key,
            "pvac_variant_id": group[0]["pvac_variant_id"],
            "gene": group[0]["gene"],
            "amino_acid_change": group[0]["amino_acid_change"],
            "hgvsp": group[0]["hgvsp"],
            "tumor_dna_depth": group[0]["tumor_dna_depth"],
            "tumor_dna_vaf": group[0]["tumor_dna_vaf"],
            "tumor_rna_depth": group[0]["tumor_rna_depth"],
            "tumor_rna_vaf": group[0]["tumor_rna_vaf"],
            "gene_expression": group[0]["gene_expression"],
            "ccf": group[0].get("ccf"),
            "clone_id": group[0].get("clone_id"),
            "clonality": group[0].get("clonality"),
            "best_hla_allele": best_binding["hla_allele"],
            "best_mt_peptide": best_binding["mt_peptide"],
            "best_mt_ic50_nm": best_binding["mt_ic50_nm"],
            "quality_2022_max_eligible": max(q22_eligible) if q22_eligible else None,
            "quality_2023_top2_mean": q23_mean,
            "quality_2023_top_peptide_1": top_two[0]["mt_peptide"] if top_two else None,
            "quality_2023_top_peptide_2": top_two[1]["mt_peptide"] if len(top_two) > 1 else None,
            "n_unique_2023_binder_peptides": len(best_per_peptide),
            "n_candidate_hla_peptide_pairs": len(group),
        }
        summaries.append(summary)
        for row in group:
            row["quality_2023_mutation_top2_mean"] = q23_mean

    summaries.sort(
        key=lambda row: (
            row["quality_2023_top2_mean"] is not None,
            row["quality_2023_top2_mean"] if row["quality_2023_top2_mean"] is not None else -math.inf,
        ),
        reverse=True,
    )
    for i, row in enumerate(summaries, 1):
        row["mutation_quality_rank"] = i
    rank_lookup = {row["mutation_id"]: row["mutation_quality_rank"] for row in summaries}
    for row in candidates:
        row["mutation_quality_rank"] = rank_lookup[row["mutation_id"]]
    return summaries


def add_joint_ranking(
    candidates: list[dict[str, Any]],
    weights: dict[str, float],
    require_pvac_filtered: bool,
) -> list[dict[str, Any]]:
    final = [
        row
        for row in candidates
        if row["eligible_2023_ic50"]
        and (not require_pvac_filtered or row["pvac_filtered_pass"])
    ]
    component_values: dict[str, list[float | None]] = {
        "q2023": [row["quality_2023_9mer"] for row in final],
        "q2022": [row["quality_2022"] for row in final],
        "binding": [
            -math.log10(row["mt_ic50_nm"])
            if row["mt_ic50_nm"] is not None and row["mt_ic50_nm"] > 0
            else None
            for row in final
        ],
        "rna_expr": [
            math.log1p(
                row["gene_expression"]
                if row["gene_expression"] is not None
                else row["transcript_expression"]
            )
            if row["gene_expression"] is not None or row["transcript_expression"] is not None
            else None
            for row in final
        ],
        "rna_vaf": [row["tumor_rna_vaf"] for row in final],
        "dna_vaf": [row["tumor_dna_vaf"] for row in final],
        "clonality": [
            row["ccf"]
            if row.get("ccf") is not None
            else (1.0 if row.get("clonality") == "clonal" else 0.0 if row.get("clonality") == "subclonal" else None)
            for row in final
        ],
    }
    percentile_components = {
        key: percentile_scores(values) for key, values in component_values.items()
    }
    for i, row in enumerate(final):
        weighted_sum = 0.0
        available_weight = 0.0
        missing: list[str] = []
        for component, weight in weights.items():
            score = percentile_components[component][i]
            row[f"component_{component}"] = score
            if score is None:
                missing.append(component)
            else:
                weighted_sum += weight * score
                available_weight += weight
        row["score_data_completeness"] = available_weight
        row["joint_score_available"] = 100.0 * weighted_sum / available_weight if available_weight else None
        row["joint_score_conservative"] = 100.0 * weighted_sum
        row["joint_score_status"] = "complete" if not missing else "provisional_missing_" + ",".join(missing)

    final.sort(
        key=lambda row: (
            row["joint_score_conservative"] if row["joint_score_conservative"] is not None else -math.inf,
            row["quality_2023_9mer"] if row["quality_2023_9mer"] is not None else -math.inf,
            -row["mt_ic50_nm"],
        ),
        reverse=True,
    )
    for rank, row in enumerate(final, 1):
        row["joint_rank"] = rank
    return final


def write_templates(output_dir: Path, summaries: list[dict[str, Any]]) -> None:
    rna_rows = [
        {
            "mutation_id": row["mutation_id"],
            "pvac_variant_id": row["pvac_variant_id"],
            "gene": row["gene"],
            "gene_expression": row["gene_expression"],
            "transcript_expression": None,
            "rna_depth": row["tumor_rna_depth"],
            "rna_vaf": row["tumor_rna_vaf"],
        }
        for row in summaries
    ]
    clone_rows = [
        {
            "mutation_id": row["mutation_id"],
            "pvac_variant_id": row["pvac_variant_id"],
            "gene": row["gene"],
            "ccf": row["ccf"],
            "clone_id": row["clone_id"],
            "clonality": row["clonality"],
        }
        for row in summaries
    ]
    write_tsv(
        output_dir / "RNA_expression_template.tsv",
        rna_rows,
        ["mutation_id", "pvac_variant_id", "gene", "gene_expression", "transcript_expression", "rna_depth", "rna_vaf"],
    )
    write_tsv(
        output_dir / "CCF_clonality_template.tsv",
        clone_rows,
        ["mutation_id", "pvac_variant_id", "gene", "ccf", "clone_id", "clonality"],
    )


def write_neoantigenediting_json(output_dir: Path, sample_id: str, candidates: list[dict[str, Any]]) -> None:
    mutations: dict[str, dict[str, Any]] = {}
    neoantigens: list[dict[str, Any]] = []
    hla_genes = sorted({row["hla_allele"].replace("HLA-", "") for row in candidates})
    for row in candidates:
        # NeoantigenEditing 2022 requires both MT and aligned WT peptides.
        # Frameshift-only rows remain in TSV outputs but are omitted here.
        if not valid_peptide(row["wt_peptide"]):
            continue
        compact_mutation_id = (
            f"{row['chromosome'].replace('chr', '')}_{row['position_1based']}_{row['reference']}_{row['variant']}"
        )
        mutations[compact_mutation_id] = {
            "id": compact_mutation_id,
            "gene": row["gene"],
            "missense": 1 if row["variant_type"].lower() == "missense" else 0,
        }
        allele_token = re.sub(r"[^A-Za-z0-9]", "", row["hla_allele"].replace("HLA-", ""))
        neo_id = (
            f"{compact_mutation_id}_{row['mutation_position_in_peptide']}_9_{allele_token}_{row['mt_peptide']}"
        )
        neoantigens.append(
            {
                "id": neo_id,
                "mutation_id": compact_mutation_id,
                "HLA_gene_id": row["hla_allele"].replace("*", ""),
                "sequence": row["mt_peptide"],
                "WT_sequence": row["wt_peptide"],
                "mutated_position": row["mutation_position_in_peptide"],
                "Kd": row["mt_ic50_nm"],
                "KdWT": row["wt_ic50_nm"],
                "R": row["q2022_R"],
                "logC": row["q2022_logC"],
                "logA": row["q2022_logA"],
                "quality": row["quality_2022"],
                "is_eligible_fitness": row["eligible_2022_ic50"],
                "quality_2023_9mer": row["quality_2023_9mer"],
            }
        )
    payload = {
        "patient": sample_id,
        "id": sample_id,
        "HLA_genes": hla_genes,
        "mutations": sorted(mutations.values(), key=lambda row: row["id"]),
        "neoantigens": neoantigens,
        "sample_trees": [],
        "note": "Peptide-level conversion only. Add PhyloWGS top trees before clone-fitness analysis.",
    }
    with (output_dir / f"{sample_id}.neoantigenediting_peptides.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pvac-tsv", required=True, type=Path, help="pVACseq MHC_I.all_epitopes.tsv")
    parser.add_argument("--aggregated-tsv", type=Path, help="Optional pVACseq all_epitopes.aggregated.tsv")
    parser.add_argument(
        "--filtered-tsv",
        type=Path,
        help="Optional pVACseq filtered.tsv; when supplied, final_candidates contains only these passing rows",
    )
    parser.add_argument("--iedb-fasta", required=True, type=Path, help="NeoantigenEditing IEDB FASTA")
    parser.add_argument("--model-json", required=True, type=Path, help="epitope_distance_model_parameters.json")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--rna-tsv", type=Path, help="Optional RNA table; use generated template columns")
    parser.add_argument("--clonality-tsv", type=Path, help="Optional CCF/clonality table; use generated template columns")
    parser.add_argument("--blastp", default="blastp")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--q2022-ic50-cutoff", type=float, default=500.0)
    parser.add_argument("--q2023-ic50-cutoff", type=float, default=4000.0)
    parser.add_argument(
        "--pvac-binding-ic50-cutoff",
        type=float,
        default=500.0,
        help="IC50 threshold used by the pVACseq run; affects pvac_binding_pass only",
    )
    parser.add_argument("--binding-percentile-cutoff", type=float, default=2.0)
    parser.add_argument("--min-tdna-depth", type=int, default=10)
    parser.add_argument("--min-tdna-vaf", type=float, default=0.05, help="Analysis flag only; does not delete rows")
    parser.add_argument("--min-normal-depth", type=int, default=5)
    parser.add_argument("--max-normal-vaf", type=float, default=0.02)
    parser.add_argument("--clonal-ccf-threshold", type=float, default=0.9)
    parser.add_argument(
        "--weights",
        default="q2023=0.25,q2022=0.15,binding=0.15,rna_expr=0.15,rna_vaf=0.10,dna_vaf=0.10,clonality=0.10",
        help="Heuristic joint-ranking weights; normalized to sum to 1",
    )
    args = parser.parse_args()

    for path in (args.pvac_tsv, args.iedb_fasta, args.model_json):
        if not path.exists():
            parser.error(f"File not found: {path}")
    if args.rna_tsv and not args.rna_tsv.exists():
        parser.error(f"RNA table not found: {args.rna_tsv}")
    if args.clonality_tsv and not args.clonality_tsv.exists():
        parser.error(f"Clonality table not found: {args.clonality_tsv}")
    if args.filtered_tsv and not args.filtered_tsv.exists():
        parser.error(f"Filtered pVACseq table not found: {args.filtered_tsv}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    weights = parse_weights(args.weights)
    pvac_rows = read_tsv(args.pvac_tsv)
    candidates = build_candidates(pvac_rows)
    for row in candidates:
        row["sample_id"] = args.sample_id

    aggregated_path = args.aggregated_tsv
    if aggregated_path is None:
        guessed = Path(str(args.pvac_tsv).replace(".all_epitopes.tsv", ".all_epitopes.aggregated.tsv"))
        aggregated_path = guessed if guessed.exists() else None
    aggregated = load_aggregated(aggregated_path)
    filtered_keys = load_filtered_candidate_keys(args.filtered_tsv)
    for row in candidates:
        agg = aggregated.get(row["pvac_variant_id"], {})
        row["pvac_tier"] = clean_text(agg.get("Tier"))
        row["pvac_evaluation"] = clean_text(agg.get("Evaluation"))
        row["pvac_aggregated_best_peptide"] = clean_text(agg.get("Best Peptide"))
        row["pvac_filtered_pass"] = (
            (
                row["pvac_variant_id"],
                row["hla_allele"],
                row["mt_peptide"],
            )
            in filtered_keys
            if args.filtered_tsv
            else row["pvac_binding_pass"] if "pvac_binding_pass" in row else False
        )

    apply_supplements(
        candidates,
        load_supplement(args.rna_tsv),
        load_supplement(args.clonality_tsv),
        args.clonal_ccf_threshold,
    )

    distance_model = EpitopeDistance(args.model_json)
    iedb_records = read_fasta(args.iedb_fasta)
    iedb_9mers = sorted(
        {
            sequence[i : i + 9]
            for sequence in iedb_records.values()
            if set(sequence) <= AA
            for i in range(max(0, len(sequence) - 8))
        }
    )
    blast_scores = blast_candidates(
        (row["mt_peptide"] for row in candidates),
        args.iedb_fasta,
        iedb_records,
        args.output_dir,
        args.blastp,
        args.makeblastdb,
    )
    add_quality_scores(candidates, distance_model, iedb_9mers, blast_scores)
    add_eligibility_and_flags(
        candidates,
        args.q2022_ic50_cutoff,
        args.q2023_ic50_cutoff,
        args.pvac_binding_ic50_cutoff,
        args.binding_percentile_cutoff,
        args.min_tdna_depth,
        args.min_tdna_vaf,
        args.min_normal_depth,
        args.max_normal_vaf,
    )
    if not args.filtered_tsv:
        for row in candidates:
            row["pvac_filtered_pass"] = row["pvac_binding_pass"]
    summaries = add_mutation_level_scores(candidates)
    final = add_joint_ranking(candidates, weights, require_pvac_filtered=bool(args.filtered_tsv))

    fields = [
        "joint_rank",
        "sample_id",
        "mutation_id",
        "pvac_variant_id",
        "chromosome",
        "position_1based",
        "start_0based",
        "stop_0based",
        "reference",
        "variant",
        "gene",
        "variant_type",
        "amino_acid_change",
        "hgvsp",
        "transcript",
        "canonical",
        "mane_select",
        "hla_allele",
        "mt_peptide",
        "wt_peptide",
        "mutation_position_in_peptide",
        "mt_ic50_nm",
        "wt_ic50_nm",
        "mt_percentile",
        "eligible_2022_ic50",
        "eligible_2023_ic50",
        "pvac_binding_pass",
        "pvac_filtered_pass",
        "q2022_R",
        "q2022_logC",
        "q2022_logA",
        "quality_2022",
        "q2023_dc_mt_wt",
        "q2023_min_dc_iedb",
        "q2023_nearest_iedb_9mer",
        "quality_2023_9mer",
        "quality_2023_mutation_top2_mean",
        "mutation_quality_rank",
        "gene_expression",
        "transcript_expression",
        "tumor_rna_depth",
        "tumor_rna_vaf",
        "tumor_dna_depth",
        "tumor_dna_vaf",
        "normal_depth",
        "normal_vaf",
        "ccf",
        "clone_id",
        "clonality",
        "clonality_source",
        "best_tree_cellular_prevalence",
        "best_tree_ccf",
        "ccf_min_top_trees",
        "ccf_max_top_trees",
        "n_trees_matched",
        "tree_assignments",
        "analysis_tdna_support_pass",
        "normal_support_pass",
        "pvac_tier",
        "pvac_evaluation",
        "component_q2023",
        "component_q2022",
        "component_binding",
        "component_rna_expr",
        "component_rna_vaf",
        "component_dna_vaf",
        "component_clonality",
        "score_data_completeness",
        "joint_score_available",
        "joint_score_conservative",
        "joint_score_status",
        "review_flags",
    ]
    all_rows = sorted(
        candidates,
        key=lambda row: (
            row["eligible_2023_ic50"],
            row.get("joint_score_conservative") if row.get("joint_score_conservative") is not None else -math.inf,
            row["quality_2023_9mer"] if row["quality_2023_9mer"] is not None else -math.inf,
            -row["mt_ic50_nm"],
        ),
        reverse=True,
    )
    write_tsv(args.output_dir / f"{args.sample_id}.all_scored.tsv", all_rows, fields)
    write_tsv(args.output_dir / f"{args.sample_id}.final_candidates.tsv", final, fields)

    summary_fields = [
        "mutation_quality_rank",
        "sample_id",
        "mutation_id",
        "pvac_variant_id",
        "gene",
        "amino_acid_change",
        "hgvsp",
        "best_hla_allele",
        "best_mt_peptide",
        "best_mt_ic50_nm",
        "quality_2022_max_eligible",
        "quality_2023_top2_mean",
        "quality_2023_top_peptide_1",
        "quality_2023_top_peptide_2",
        "n_unique_2023_binder_peptides",
        "n_candidate_hla_peptide_pairs",
        "gene_expression",
        "tumor_rna_depth",
        "tumor_rna_vaf",
        "tumor_dna_depth",
        "tumor_dna_vaf",
        "ccf",
        "clone_id",
        "clonality",
    ]
    write_tsv(args.output_dir / f"{args.sample_id}.mutation_summary.tsv", summaries, summary_fields)
    write_templates(args.output_dir, summaries)
    write_neoantigenediting_json(args.output_dir, args.sample_id, candidates)

    manifest = {
        "sample_id": args.sample_id,
        "hla_alleles": sorted({row["hla_allele"] for row in candidates}),
        "pvac_input": str(args.pvac_tsv),
        "aggregated_input": str(aggregated_path) if aggregated_path else None,
        "filtered_input": str(args.filtered_tsv) if args.filtered_tsv else None,
        "iedb_fasta": str(args.iedb_fasta),
        "epitope_distance_model": str(args.model_json),
        "source_rows": len(pvac_rows),
        "deduplicated_hla_peptide_candidates": len(candidates),
        "final_2023_binder_candidates": len(final),
        "mutations": len(summaries),
        "rna_input": str(args.rna_tsv) if args.rna_tsv else None,
        "clonality_input": str(args.clonality_tsv) if args.clonality_tsv else None,
        "weights": weights,
        "thresholds": {
            "q2022_ic50_nm": args.q2022_ic50_cutoff,
            "q2023_ic50_nm": args.q2023_ic50_cutoff,
            "pvac_binding_ic50_nm": args.pvac_binding_ic50_cutoff,
            "binding_percentile": args.binding_percentile_cutoff,
            "analysis_min_tdna_depth": args.min_tdna_depth,
            "analysis_min_tdna_vaf": args.min_tdna_vaf,
            "min_normal_depth": args.min_normal_depth,
            "max_normal_vaf": args.max_normal_vaf,
            "clonal_ccf": args.clonal_ccf_threshold,
        },
        "ranking_note": "Joint ranking is a transparent heuristic, not a published clinical classifier. Missing components are not imputed; completeness is reported.",
    }
    with (args.output_dir / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)

    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"External command failed: {exc}", file=sys.stderr)
        if exc.stderr:
            print(exc.stderr, file=sys.stderr)
        raise SystemExit(2)
