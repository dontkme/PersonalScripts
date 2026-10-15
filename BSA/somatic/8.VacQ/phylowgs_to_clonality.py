#!/usr/bin/env python3
"""Map pVACseq variants to PhyloWGS top-tree clones and estimate tumor CCF.

PhyloWGS cellular prevalence is a fraction of all cells in the sequenced
sample.  For a tumor/normal mixture, this script additionally reports a
purity-normalized cancer-cell fraction (CCF = cellular prevalence / purity,
capped at 1).  Clone labels are tree-specific, so assignments retain tree IDs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any


MISSING = {"", "NA", "N/A", "NONE", "NULL", "."}


def clean(value: Any) -> str:
    if value is None or str(value).strip().upper() in MISSING:
        return ""
    return str(value).strip()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = {}
            for field in fields:
                value = row.get(field)
                if value is None:
                    output[field] = "NA"
                elif isinstance(value, float):
                    output[field] = f"{value:.8g}"
                else:
                    output[field] = value
            writer.writerow(output)


def normalize_chromosome(value: str) -> str:
    value = clean(value)
    return value[3:] if value.lower().startswith("chr") else value


def load_ssms(path: Path) -> tuple[dict[str, tuple[str, int]], dict[tuple[str, int], list[str]]]:
    by_id: dict[str, tuple[str, int]] = {}
    by_coordinate: dict[tuple[str, int], list[str]] = defaultdict(list)
    for row in read_tsv(path):
        token = clean(row.get("gene"))
        chrom, separator, position_text = token.rpartition("_")
        if not separator:
            raise ValueError(f"Cannot parse PhyloWGS SSM coordinate: {token}")
        coordinate = (normalize_chromosome(chrom), int(position_text))
        ssm_id = clean(row.get("id"))
        by_id[ssm_id] = coordinate
        by_coordinate[coordinate].append(ssm_id)
    return by_id, by_coordinate


def load_tree_order(path: Path | None, summary: dict[str, Any]) -> list[str]:
    if path is not None:
        rows = read_tsv(path)
        return [clean(row.get("tree_id")) for row in rows if clean(row.get("tree_id"))]
    return list(summary["trees"].keys())


def mutation_id(row: dict[str, str]) -> tuple[str, str, int]:
    chrom_raw = clean(row.get("Chromosome"))
    chrom = normalize_chromosome(chrom_raw)
    position = int(float(clean(row.get("Start")))) + 1
    ref = clean(row.get("Reference"))
    alt = clean(row.get("Variant"))
    return f"{chrom_raw}:{position}:{ref}>{alt}", chrom, position


def pvac_variant_id(row: dict[str, str]) -> str:
    return "-".join(
        clean(row.get(field))
        for field in ("Chromosome", "Start", "Stop", "Reference", "Variant")
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--ssm-data", required=True, type=Path)
    parser.add_argument("--summary-gz", required=True, type=Path)
    parser.add_argument("--mutass-zip", required=True, type=Path)
    parser.add_argument("--pvac-tsv", required=True, type=Path)
    parser.add_argument("--tree-ids", type=Path, help="Optional ranked tree_id TSV")
    parser.add_argument("--tumor-purity", required=True, type=float)
    parser.add_argument("--clonal-ccf-threshold", type=float, default=0.9)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if not 0 < args.tumor_purity <= 1:
        parser.error("--tumor-purity must be in (0, 1]")

    _, ssms_by_coordinate = load_ssms(args.ssm_data)
    with gzip.open(args.summary_gz, "rt", encoding="utf-8") as handle:
        summary = json.load(handle)
    tree_order = load_tree_order(args.tree_ids, summary)
    if not tree_order:
        raise ValueError("No tree IDs found")

    tree_llh = {tree_id: float(summary["trees"][tree_id]["llh"]) for tree_id in tree_order}
    max_llh = max(tree_llh.values())
    raw_weights = {tree_id: math.exp(tree_llh[tree_id] - max_llh) for tree_id in tree_order}
    weight_sum = sum(raw_weights.values())
    tree_weights = {tree_id: value / weight_sum for tree_id, value in raw_weights.items()}

    assignments: dict[str, dict[str, str]] = {}
    with zipfile.ZipFile(args.mutass_zip) as archive:
        for tree_id in tree_order:
            payload = json.loads(archive.read(f"{tree_id}.json"))
            assignments[tree_id] = {}
            for clone_id, members in payload["mut_assignments"].items():
                for ssm_id in members.get("ssms", []):
                    assignments[tree_id][ssm_id] = clone_id

    unique_variants: dict[str, dict[str, str]] = {}
    for row in read_tsv(args.pvac_tsv):
        variant_id = pvac_variant_id(row)
        previous = unique_variants.get(variant_id)
        if previous is None or (not clean(previous.get("Gene Name")) and clean(row.get("Gene Name"))):
            unique_variants[variant_id] = row

    output_rows: list[dict[str, Any]] = []
    for variant_id, pvac_row in unique_variants.items():
        mut_id, chrom, position = mutation_id(pvac_row)
        ssm_ids = ssms_by_coordinate.get((chrom, position), [])
        row: dict[str, Any] = {
            "mutation_id": mut_id,
            "pvac_variant_id": variant_id,
            "gene": clean(pvac_row.get("Gene Name")),
            "chromosome": clean(pvac_row.get("Chromosome")),
            "position_1based": position,
            "reference": clean(pvac_row.get("Reference")),
            "variant": clean(pvac_row.get("Variant")),
            "ssm_id": ",".join(ssm_ids),
        }
        per_tree: list[dict[str, Any]] = []
        for tree_id in tree_order:
            clone_ids = sorted({assignments[tree_id][ssm] for ssm in ssm_ids if ssm in assignments[tree_id]})
            if len(clone_ids) != 1:
                continue
            clone_id = clone_ids[0]
            prevalence = float(summary["trees"][tree_id]["populations"][clone_id]["cellular_prevalence"][0])
            ccf = min(1.0, prevalence / args.tumor_purity)
            per_tree.append(
                {
                    "tree_id": tree_id,
                    "clone_id": clone_id,
                    "prevalence": prevalence,
                    "ccf": ccf,
                    "weight": tree_weights[tree_id],
                }
            )

        if per_tree:
            best = next(item for item in per_tree if item["tree_id"] == tree_order[0])
            present_weight = sum(item["weight"] for item in per_tree)
            weighted_cp = sum(item["prevalence"] * item["weight"] for item in per_tree) / present_weight
            weighted_ccf = sum(item["ccf"] * item["weight"] for item in per_tree) / present_weight
            simple_ccfs = [item["ccf"] for item in per_tree]
            row.update(
                {
                    "match_status": "matched",
                    "clone_id": f"{best['tree_id']}:{best['clone_id']}",
                    "best_tree_id": best["tree_id"],
                    "best_tree_clone_id": best["clone_id"],
                    "best_tree_cellular_prevalence": best["prevalence"],
                    "best_tree_ccf": best["ccf"],
                    "cellular_prevalence_weighted_top_trees": weighted_cp,
                    "ccf": weighted_ccf,
                    "ccf_mean_top_trees": sum(simple_ccfs) / len(simple_ccfs),
                    "ccf_min_top_trees": min(simple_ccfs),
                    "ccf_max_top_trees": max(simple_ccfs),
                    "n_trees_matched": len(per_tree),
                    "clonality": "clonal" if weighted_ccf >= args.clonal_ccf_threshold else "subclonal",
                    "tree_assignments": ";".join(
                        f"{item['tree_id']}:clone{item['clone_id']}:CP={item['prevalence']:.6g}:CCF={item['ccf']:.6g}"
                        for item in per_tree
                    ),
                }
            )
        else:
            row.update(
                {
                    "match_status": "not_in_phylowgs_ssm",
                    "clone_id": "",
                    "best_tree_id": "",
                    "best_tree_clone_id": "",
                    "best_tree_cellular_prevalence": None,
                    "best_tree_ccf": None,
                    "cellular_prevalence_weighted_top_trees": None,
                    "ccf": None,
                    "ccf_mean_top_trees": None,
                    "ccf_min_top_trees": None,
                    "ccf_max_top_trees": None,
                    "n_trees_matched": 0,
                    "clonality": "",
                    "tree_assignments": "",
                }
            )
        output_rows.append(row)

    fields = [
        "mutation_id",
        "pvac_variant_id",
        "gene",
        "chromosome",
        "position_1based",
        "reference",
        "variant",
        "ssm_id",
        "match_status",
        "clone_id",
        "best_tree_id",
        "best_tree_clone_id",
        "best_tree_cellular_prevalence",
        "best_tree_ccf",
        "cellular_prevalence_weighted_top_trees",
        "ccf",
        "ccf_mean_top_trees",
        "ccf_min_top_trees",
        "ccf_max_top_trees",
        "n_trees_matched",
        "clonality",
        "tree_assignments",
    ]
    output_rows.sort(key=lambda row: (normalize_chromosome(row["chromosome"]), row["position_1based"]))
    write_tsv(args.output, output_rows, fields)

    manifest = {
        "pVAC_variants": len(output_rows),
        "matched_to_PhyloWGS_SSM": sum(row["match_status"] == "matched" for row in output_rows),
        "not_in_PhyloWGS_SSM": sum(row["match_status"] != "matched" for row in output_rows),
        "tree_order": tree_order,
        "relative_top_tree_weights": tree_weights,
        "tumor_purity": args.tumor_purity,
        "clonal_ccf_threshold": args.clonal_ccf_threshold,
    }
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
