#!/usr/bin/env python3
"""Aggregate OE_WGS pipeline outputs into per-patient and cohort reports.

Reads existing pipeline outputs only -- never re-runs any analysis.

    python3 make_reports.py [--base-dir .] [--patients P ...] [--outdir reports]
                            [--genes-dir CancerGens] [--count-method auto|bcftools|gzip]
                            [--verify] [--strict]

Writes reports/<P>.report.md, reports/cohort_summary.tsv, reports/SUMMARY.md
and reports/tables/*.tsv.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FLAGSTAT_RE = re.compile(
    r"^(?P<n>\d+) \+ (?P<nf>\d+) (?P<label>.+?)"
    r"(?: \((?P<pct>[\d.]+)% : N/A\))?$"
)

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

# Effect sets matching the SnpSift filters used to build each TSV. A site only
# enters a TSV because SOME transcript matched, so per-transcript aggregation
# must be limited to the matching transcripts -- otherwise a gene whose only
# annotation is e.g. downstream_gene_variant gets reported as a hit.
PROTEIN_DAMAGING_EFFECTS = {
    "missense_variant", "frameshift_variant", "stop_gained",
    "splice_acceptor_variant", "splice_donor_variant",
    "start_lost", "stop_lost",
}
MISSENSE_FRAMESHIFT_EFFECTS = {"missense_variant", "frameshift_variant"}


def effect_matches(effect: str, targets: set[str]) -> bool:
    """True if any '&'-joined component of a compound effect is targeted."""
    return any(part.strip() in targets for part in effect.split("&"))

# (side, subdir, sample-prefix template)
SIDES = [
    ("tumor", "WGS", "PC_{p}.Tumor"),
    ("blood", "WGS_blood", "{p}.Blood"),
]

CHAIN_LABELS = {
    "mutect2": "Mutect2-only",
    "both": "Mutect2 + Strelka2 both PASS",
}

# Per-chain filename stems, relative to 5.FilterSNV/
CHAIN_STEMS = {
    "mutect2": "mutect2.filtered.PASS.norm",
    "both": "mutect2_strelka.both.PASS",
}

ANN_HEADER = [
    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
    "ANN[*].GENE", "ANN[*].EFFECT", "ANN[*].IMPACT", "ANN[*].FEATUREID",
    "ANN[*].HGVS_C", "ANN[*].HGVS_P", "ANN[*].RANK", "TLOD",
]

NULL_MARK = "n/a"
# Top-level dirs that hold shared assets rather than a patient
SKIP_DIRS = {"WGS", "WGS_blood", "CancerGens", "reference", "reports",
             "6.1phyloWGS", "8.VacQ", "NeoantigenEditing-main"}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def find_one(root: Path, pattern: str) -> Path | None:
    if not root.is_dir():
        return None
    hits = sorted(root.glob(pattern))
    return hits[0] if hits else None


def read_tsv(path: Path) -> list[dict]:
    """Read a TSV into a list of dicts (header row required)."""
    if path is None or not path.exists():
        return []
    with open(path, errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def safe_int(value) -> int | None:
    f = safe_float(value)
    return int(f) if f is not None else None


def fmt_int(n) -> str:
    return f"{n:,}" if isinstance(n, int) else NULL_MARK


def fmt_float(x, digits: int = 4) -> str:
    if x is None:
        return NULL_MARK
    if isinstance(x, float) and (abs(x) < 1e-3 and x != 0):
        return f"{x:.3e}"
    return f"{x:.{digits}f}"


def fmt_pct(x, digits: int = 2) -> str:
    return NULL_MARK if x is None else f"{x * 100:.{digits}f}%"


def median(values) -> float | None:
    vals = [v for v in values if v is not None]
    return statistics.median(vals) if vals else None


def md_table(headers: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# VCF counting (bcftools preferred, gzip fallback)
# ---------------------------------------------------------------------------

def vcf_stats(path: Path | None, method: str = "auto",
              bcftools: str | None = None) -> tuple[int | None, Counter]:
    """Return (n_records, Counter(filter_column)) for a VCF, or (None, {})."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return None, Counter()

    lines: list[str] | None = None
    if method in ("auto", "bcftools") and bcftools:
        try:
            proc = subprocess.run([bcftools, "view", "-H", str(path)],
                                  capture_output=True, check=False)
            if proc.returncode == 0:
                lines = proc.stdout.decode("utf-8", "replace").splitlines()
        except OSError:
            lines = None
        if lines is None and method == "bcftools":
            return None, Counter()

    if lines is None:
        try:
            with gzip.open(path, "rt", errors="replace") as fh:
                lines = [ln for ln in fh if not ln.startswith("#")]
        except (gzip.BadGzipFile, OSError, EOFError):
            try:
                with open(path, errors="replace") as fh:
                    lines = [ln for ln in fh if not ln.startswith("#")]
            except OSError:
                return None, Counter()

    tally: Counter = Counter()
    count = 0
    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue
        count += 1
        fields = line.split("\t")
        if len(fields) > 6:
            tally[fields[6]] += 1
    return count, tally


def count_vcf_records(path: Path | None, method: str = "auto",
                      bcftools: str | None = None) -> int | None:
    return vcf_stats(path, method, bcftools)[0]


# ---------------------------------------------------------------------------
# Stage parsers
# ---------------------------------------------------------------------------

def parse_flagstat(path: Path | None) -> dict:
    """Parse samtools flagstat. Keys on label text, never line index."""
    if path is None or not path.exists():
        return {}
    out: dict = {}
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                m = FLAGSTAT_RE.match(line.strip())
                if not m:
                    continue
                # The non-percent parenthetical on the "in total" line
                # ("(QC-passed reads + QC-failed reads)") lands in the label;
                # strip any trailing bracket group to get the bare label.
                label = re.sub(r"\s*\([^)]*\)\s*$", "", m.group("label")).strip()
                pct = safe_float(m.group("pct"))
                # First occurrence wins; "mapped" must not be clobbered by
                # "primary mapped" etc.
                if label not in out:
                    out[label] = {"count": int(m.group("n")), "pct": pct}
    except OSError:
        return {}
    return out


def parse_picard_metrics(path: Path | None) -> dict:
    """Parse a Picard MarkDuplicates metrics file (tolerates histogram block)."""
    if path is None or not path.exists():
        return {}
    header = None
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                if line.startswith("## HISTOGRAM"):
                    break
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if header is None:
                    if "PERCENT_DUPLICATION" in fields:
                        header = fields
                    continue
                return dict(zip(header, fields))
    except OSError:
        return {}
    return {}


def parse_fastp_json(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {}
    try:
        with open(path, errors="replace") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}

    def sect(name, key):
        return data.get(name, {}).get(key)

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    return {
        "raw_reads": before.get("total_reads"),
        "clean_reads": after.get("total_reads"),
        "raw_bases": before.get("total_bases"),
        "clean_bases": after.get("total_bases"),
        "q30_before": before.get("q30_rate"),
        "q30_after": after.get("q30_rate"),
        "gc_before": before.get("gc_content"),
        "gc_after": after.get("gc_content"),
        "passed_filter_reads": sect("filtering_result", "passed_filter_reads"),
        "low_quality_reads": sect("filtering_result", "low_quality_reads"),
        "too_short_reads": sect("filtering_result", "too_short_reads"),
        "adapter_dimer_reads": sect("filtering_result", "adapter_dimer_reads"),
        "dup_rate": sect("duplication", "rate"),
        "adapter_trimmed_reads": sect("adapter_cutting", "adapter_trimmed_reads"),
        "insert_peak": sect("insert_size", "peak"),
    }


def parse_filtering_stats(path: Path | None) -> dict:
    """Parse GATK --filtering-stats output."""
    out: dict = {"meta": {}, "filters": []}
    if path is None or not path.exists():
        return out
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#<METADATA>"):
                    body = line[len("#<METADATA>"):]
                    if "=" in body:
                        key, val = body.split("=", 1)
                        key = key.strip()
                        if key in ("threshold", "fdr", "sensitivity"):
                            out["meta"][key] = safe_float(val)
                    continue
                if line.startswith("filter\t"):
                    continue
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) >= 5:
                    out["filters"].append({
                        "filter": fields[0],
                        "FP": safe_float(fields[1]),
                        "FDR": safe_float(fields[2]),
                        "FN": safe_float(fields[3]),
                        "FNR": safe_float(fields[4]),
                    })
    except OSError:
        return out
    return out


def parse_annotation_tsv(path: Path | None,
                         targets: set[str] | None = None) -> dict:
    """Parse a snpEff extractFields TSV.

    GENE / EFFECT / IMPACT are comma-separated parallel lists (one entry per
    transcript) within a single cell. Index i of each list belongs together.

    `targets` limits aggregation to transcripts whose effect matches the
    SnpSift filter that produced the TSV, so genes carried along only as
    intronic/downstream neighbours are not counted as hits.
    """
    result = {
        "path": path, "n_sites": 0, "impacts": Counter(), "effects": Counter(),
        "effect_any": Counter(), "tlods": [], "genes": {}, "readable": False,
        "header_ok": False, "n_transcripts_matched": 0,
    }
    if path is None or not path.exists():
        return result

    with open(path, errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return result
        result["header_ok"] = reader.fieldnames == ANN_HEADER
        gcol, ecol, icol = "ANN[*].GENE", "ANN[*].EFFECT", "ANN[*].IMPACT"
        if gcol not in reader.fieldnames:
            return result
        # A parseable header means the filter ran; 0 data rows is a valid
        # "no variants passed" result, distinct from a missing file.
        result["readable"] = True
        for row in reader:
            result["n_sites"] += 1
            site = (row.get("CHROM"), row.get("POS"), row.get("REF"), row.get("ALT"))
            tlod = safe_float(row.get("TLOD"))
            if tlod is not None:
                result["tlods"].append(tlod)

            genes = (row.get(gcol) or "").split(",")
            effs = (row.get(ecol) or "").split(",")
            imps = (row.get(icol) or "").split(",")

            seen_genes: set[str] = set()
            site_effects: set[str] = set()
            for i in range(min(len(genes), len(effs), len(imps))):
                gene = genes[i].strip()
                eff = effs[i].strip()
                imp = imps[i].strip()
                if targets is not None and not effect_matches(eff, targets):
                    continue
                result["n_transcripts_matched"] += 1
                result["impacts"][imp] += 1
                for part in eff.split("&"):
                    part = part.strip()
                    if part:
                        result["effect_any"][part] += 1
                        site_effects.add(part)
                if not gene or gene in ("NA", "."):
                    continue
                if gene in seen_genes:
                    continue
                seen_genes.add(gene)

                entry = result["genes"].setdefault(gene, {
                    "sites": set(), "impacts": Counter(), "effects": Counter(),
                    "tlod_max": None, "worst": None, "high_sites": set(),
                })
                entry["sites"].add(site)
                entry["impacts"][imp] += 1
                entry["effects"][eff] += 1
                if imp == "HIGH":
                    entry["high_sites"].add(site)
                rank = IMPACT_RANK.get(imp, 9)
                if entry["worst"] is None or rank < entry["worst"][0]:
                    entry["worst"] = (rank, imp, eff)
                if tlod is not None and (entry["tlod_max"] is None or tlod > entry["tlod_max"]):
                    entry["tlod_max"] = tlod

            for eff in site_effects:
                result["effects"][eff] += 1

    return result


def resolve_strelka_vcf(patient_dir: Path, kind: str) -> Path | None:
    """Locate 4_1.St_SNV somatic.{kind}.vcf.gz.

    Nested run dirs differ per patient (E88354 uses 4_1.St_SNV/4.1St/...),
    so glob recursively. Skip any *_1K* run so the gnomAD variant doesn't win.
    """
    root = patient_dir / "4_1.St_SNV"
    if not root.is_dir():
        return None
    hits = []
    for path in root.rglob(f"somatic.{kind}.vcf.gz"):
        if any(part.lower().endswith("_1k") or "1k" == part.lower()
               for part in path.parts):
            continue
        hits.append(path)
    if not hits:
        return None
    hits.sort(key=lambda p: (len(p.parts), -p.stat().st_mtime))
    return hits[0]


def resolve_hla_bestguess(patient_dir: Path, patient: str) -> tuple[Path | None, list[Path]]:
    """Find hla/R1_bestguess_G.txt, handling BE0728's nested layout."""
    root = patient_dir / "7.HLA_type"
    if not root.is_dir():
        return None, []
    hits = sorted(root.glob("**/hla/R1_bestguess_G.txt"))
    if not hits:
        return None, []
    if len(hits) == 1:
        return hits[0], hits
    preferred = [p for p in hits
                 if "IGV" not in str(p) and patient.lower() in str(p).lower()]
    return (preferred[0] if preferred else hits[0]), hits


def parse_hla(path: Path | None) -> dict:
    """Pivot R1_bestguess_G.txt (two rows per locus) into per-locus records."""
    out = {"loci": [], "n_rows": 0, "n_loci": 0, "n_perfect_g": 0, "n_ambiguous": 0}
    rows = read_tsv(path)
    if not rows:
        return out
    out["n_rows"] = len(rows)

    by_locus: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_locus[(row.get("Locus") or "").strip()].append(row)

    for locus in sorted(by_locus):
        entries = sorted(by_locus[locus],
                         key=lambda r: safe_int(r.get("Chromosome")) or 0)
        alleles, q1s, covs, perfect = [], [], [], []
        ambiguous = False
        for entry in entries:
            raw = (entry.get("Allele") or "").strip()
            allele_list = [a.strip() for a in raw.split(";") if a.strip()]
            if len(allele_list) > 1:
                ambiguous = True
            alleles.append(allele_list)
            q1s.append(safe_float(entry.get("Q1")))
            covs.append(safe_float(entry.get("AverageCoverage")))
            perfect.append(safe_int(entry.get("perfectG")))
        n_perfect = sum(1 for p in perfect if p == 1)
        out["n_perfect_g"] += n_perfect
        if ambiguous:
            out["n_ambiguous"] += 1
        out["loci"].append({
            "locus": locus,
            "alleles": alleles,
            "q1": q1s,
            "coverage": covs,
            "perfect": perfect,
            "n_perfect": n_perfect,
            "ambiguous": ambiguous,
        })
    out["n_loci"] = len(out["loci"])
    return out


# ---------------------------------------------------------------------------
# Cancer gene lists
# ---------------------------------------------------------------------------

def load_cancer_genes(genes_dir: Path) -> tuple[dict, dict]:
    """Merge COSMIC / IntOGen / OncoKB lists into {SYMBOL: {...}}."""
    merged: dict[str, dict] = {}
    stats = {"sources": {}, "union": 0, "all_three": 0, "files": []}

    def entry(symbol: str) -> dict:
        key = symbol.strip().upper()
        return merged.setdefault(key, {"sources": set(), "intogen_roles": set(),
                                       "intogen_cohorts": 0})

    def load(filename: str, source: str, symbol_col: str, extra=None):
        path = genes_dir / filename
        if not path.exists():
            return
        stats["files"].append(str(path))
        symbols = set()
        for row in read_tsv(path):
            symbol = (row.get(symbol_col) or "").strip()
            if not symbol:
                continue
            symbols.add(symbol.upper())
            item = entry(symbol)
            item["sources"].add(source)
            if extra:
                extra(item, row)
        stats["sources"][source] = len(symbols)

    load("cancerGeneList.tsv", "OncoKB", "Hugo Symbol",
         lambda it, r: it.update({
             "oncokb_gene_type": (r.get("Gene Type") or "").strip(),
             "oncokb_annotated": (r.get("OncoKB Annotated") or "").strip(),
         }))
    load("Compendium_Cancer_Genes.tsv", "IntOGen", "SYMBOL",
         lambda it, r: (it["intogen_roles"].add((r.get("ROLE") or "").strip()),
                        it.update(intogen_cohorts=it["intogen_cohorts"] + 1)))
    load("Cosmic_CancerGeneCensus_v104_GRCh38.tsv", "COSMIC", "GENE_SYMBOL",
         lambda it, r: it.update({
             "cosmic_tier": safe_int(r.get("TIER")),
             "cosmic_role": (r.get("ROLE_IN_CANCER") or "").strip(),
             "cosmic_somatic": (r.get("SOMATIC") or "").strip(),
         }))

    for item in merged.values():
        item["sources"] = sorted(item["sources"])
        item["tier"] = len(item["sources"])
        item["intogen_roles"] = sorted(r for r in item["intogen_roles"] if r)

    stats["union"] = len(merged)
    stats["all_three"] = sum(1 for i in merged.values() if i["tier"] == 3)
    return merged, stats


def cancer_flag(gene: str, genes: dict) -> dict | None:
    return genes.get(gene.strip().upper())


def cancer_cell(info: dict | None) -> str:
    if not info:
        return "no"
    srcs = ",".join(info["sources"])
    detail = []
    if info.get("cosmic_tier"):
        detail.append(f"COSMIC T{info['cosmic_tier']}")
    if info.get("oncokb_gene_type"):
        detail.append(f"OncoKB {info['oncokb_gene_type']}")
    if info.get("intogen_roles"):
        detail.append("IntOGen " + "/".join(info["intogen_roles"]))
    suffix = f" ({'; '.join(detail)})" if detail else ""
    return f"**yes** [{info['tier']}/3: {srcs}]{suffix}"


# ---------------------------------------------------------------------------
# Per-patient extraction
# ---------------------------------------------------------------------------

def extract_patient(patient: str, base: Path, genes: dict,
                    method: str, bcftools: str | None) -> dict:
    pdir = base / patient
    data: dict = {"patient": patient, "dir": pdir, "notes": [],
                  "sides": {}, "calling": {}, "filtering": {}, "chains": {},
                  "hla": {}, "genes": [], "status": {}}

    # ---- QC + alignment + BQSR, per side ----
    for side, subdir, prefix_tpl in SIDES:
        prefix = prefix_tpl.format(p=patient)
        sroot = pdir / subdir
        qc = parse_fastp_json(find_one(sroot / "1.QC", f"{prefix}.fastp_report.json"))
        if qc.get("raw_reads") and qc.get("clean_reads"):
            qc["pass_rate"] = qc["clean_reads"] / qc["raw_reads"]
        raw_fs = parse_flagstat(find_one(sroot / "2.BWA2_res", f"{prefix}.flagstat"))
        mark_fs = parse_flagstat(find_one(sroot / "2.BWA2_res", f"{prefix}.markdup.flagstat"))
        rec_fs = parse_flagstat(find_one(sroot / "3.BQSR", f"{prefix}.markdup.recal.bam.flagstat"))
        picard = parse_picard_metrics(find_one(sroot / "2.BWA2_res", f"{prefix}.rmdup.metrics.txt"))
        dup_pct = safe_float(picard.get("PERCENT_DUPLICATION"))
        if dup_pct is None and mark_fs.get("duplicates") and mark_fs.get("in total"):
            total = mark_fs["in total"]["count"]
            dup_pct = mark_fs["duplicates"]["count"] / total if total else None
        data["sides"][side] = {
            "prefix": prefix, "qc": qc, "raw_flagstat": raw_fs,
            "markdup_flagstat": mark_fs, "recal_flagstat": rec_fs,
            "picard": picard, "dup_pct": dup_pct,
        }

    # ---- Contamination ----
    contam_rows = read_tsv(find_one(pdir / "3.2CalCont", f"{patient}.contamination.table"))
    data["contamination"] = {
        "contamination": safe_float(contam_rows[0].get("contamination")) if contam_rows else None,
        "error": safe_float(contam_rows[0].get("error")) if contam_rows else None,
    }

    # ---- Calling ----
    unfiltered = find_one(pdir / "4.SNV", f"{patient}.recal.somatic.SNV.unfiltered.vcf.gz")
    n_unfiltered = count_vcf_records(unfiltered, method, bcftools)
    strelka = {}
    for kind in ("snvs", "indels"):
        path = resolve_strelka_vcf(pdir, kind)
        n, tally = vcf_stats(path, method, bcftools)
        strelka[kind] = {"path": path, "n": n, "filters": tally}
    data["calling"] = {"unfiltered": n_unfiltered, "unfiltered_path": unfiltered,
                       "strelka": strelka}

    # ---- Filtering ----
    fs_path = find_one(pdir / "5.FilterSNV", f"{patient}.filtering.stats")
    mut_pass = find_one(pdir / "5.FilterSNV", f"{patient}.mutect2.filtered.PASS.norm.vcf.gz")
    both_pass = find_one(pdir / "5.FilterSNV", f"{patient}.mutect2_strelka.both.PASS.vcf.gz")
    mut_tagged = find_one(pdir / "5.FilterSNV",
                          f"{patient}.mutect2.filtered.PASS.norm.strelka.vcf.gz")
    n_mut = count_vcf_records(mut_pass, method, bcftools)
    n_both = count_vcf_records(both_pass, method, bcftools)
    n_tagged, tag_tally = vcf_stats(mut_tagged, method, bcftools)
    data["filtering"] = {
        "stats": parse_filtering_stats(fs_path), "stats_path": fs_path,
        "mutect2_pass": n_mut, "mutect2_pass_path": mut_pass,
        "both_pass": n_both, "both_pass_path": both_pass,
        "tagged": n_tagged, "tag_tally": tag_tally,
    }

    # ---- Annotation chains ----
    for chain, stem in CHAIN_STEMS.items():
        mf = find_one(pdir / "5.FilterSNV", f"{patient}.{stem}.snpeff.missense_frameshift.tsv")
        pd = find_one(pdir / "5.FilterSNV", f"{patient}.{stem}.snpeff.protein_damaging.tsv")
        vcf = find_one(pdir / "5.FilterSNV", f"{patient}.{stem}.snpeff.vcf")
        data["chains"][chain] = {
            "stem": stem,
            "missense_frameshift": parse_annotation_tsv(mf, MISSENSE_FRAMESHIFT_EFFECTS),
            "protein_damaging": parse_annotation_tsv(pd, PROTEIN_DAMAGING_EFFECTS),
            "vcf": vcf,
            "vcf_size": vcf.stat().st_size if vcf and vcf.exists() else 0,
        }

    # ---- Gene table from the protein_damaging chain (mutect2-only) ----
    primary = data["chains"]["mutect2"]["protein_damaging"]
    gene_rows = []
    for gene, entry in primary["genes"].items():
        info = cancer_flag(gene, genes)
        gene_rows.append({
            "gene": gene,
            "n_sites": len(entry["sites"]),
            "high_sites": len(entry["high_sites"]),
            "worst_impact": entry["worst"][1] if entry["worst"] else "",
            "worst_consequence": entry["worst"][2] if entry["worst"] else "",
            "tlod_max": entry["tlod_max"],
            "cancer": info,
        })
    gene_rows.sort(key=lambda r: (-(r["cancer"]["tier"] if r["cancer"] else 0),
                                  -r["n_sites"],
                                  IMPACT_RANK.get(r["worst_impact"], 9),
                                  r["gene"]))
    data["genes"] = gene_rows
    data["gene_source"] = "protein_damaging (Mutect2-only)"

    # ---- HLA ----
    hla_path, hla_hits = resolve_hla_bestguess(pdir, patient)
    if len(hla_hits) > 1:
        data["notes"].append(f"Multiple R1_bestguess_G.txt candidates: "
                             f"{', '.join(str(h) for h in hla_hits)}")
    data["hla"] = parse_hla(hla_path)
    data["hla"]["path"] = hla_path

    # ---- Legacy fallback + completeness ----
    data["status"] = check_completeness(data)
    return data


def check_completeness(data: dict) -> dict:
    """Derive per-stage status from the data itself (never a hardcoded list)."""
    status = {}

    def sides_ok(key):
        seen = [s["qc"].get(key) for s in data["sides"].values()]
        return all(v for v in seen) if seen else False

    status["qc"] = ("ok", "") if sides_ok("clean_reads") else (
        "missing", "fastp report missing or empty")
    status["alignment"] = ("ok", "") if all(
        s["markdup_flagstat"].get("mapped") for s in data["sides"].values()
    ) else ("missing", "markdup flagstat missing")
    status["bqsr"] = ("ok", "") if all(
        s["recal_flagstat"].get("mapped") for s in data["sides"].values()
    ) else ("missing", "recal flagstat missing")
    status["contamination"] = ("ok", "") if data["contamination"]["contamination"] is not None \
        else ("missing", "contamination table missing")

    calling = data["calling"]
    if calling["unfiltered"]:
        status["calling"] = ("ok", "")
    else:
        status["calling"] = ("missing", "Mutect2 unfiltered VCF missing or empty")

    filt = data["filtering"]
    if filt["mutect2_pass"] and filt["stats"]["filters"]:
        status["filtering"] = ("ok", "")
    elif filt["mutect2_pass"]:
        status["filtering"] = ("incomplete", "filtering.stats missing")
    else:
        status["filtering"] = ("incomplete", "Mutect2 PASS VCF has 0 records")

    if filt["mutect2_pass"] and not filt["both_pass"]:
        status["filtering"] = ("incomplete",
                               "Strelka2 merge produced 0 records")

    chain_status = {}
    for chain, info in data["chains"].items():
        mf, pd = info["missense_frameshift"], info["protein_damaging"]
        if mf["readable"] or pd["readable"]:
            # A readable header means the chain ran; 0 rows is a valid result.
            chain_status[chain] = ("ok", "")
        elif info["vcf_size"] == 0:
            chain_status[chain] = ("missing", "snpEff outputs absent")
        else:
            chain_status[chain] = ("missing", "snpEff TSVs unreadable")
    status["annotation"] = chain_status.get("mutect2", ("missing", ""))

    # Stale: step-5 predates the upstream Mutect2 VCF, or is empty upstream-populated
    unfiltered_path = data["calling"]["unfiltered_path"]
    pass_path = data["filtering"]["mutect2_pass_path"]
    if unfiltered_path and unfiltered_path.exists():
        up_mtime = unfiltered_path.stat().st_mtime
        candidates = [p for p in (pass_path, data["chains"]["mutect2"]["vcf"]) if p]
        older = [p for p in candidates if p.exists() and p.stat().st_mtime < up_mtime]
        if older:
            status["filtering"] = ("stale", "step-5 output older than Mutect2 VCF")
    if data["filtering"]["mutect2_pass"] == 0 and data["calling"]["unfiltered"]:
        status["filtering"] = ("incomplete",
                               "step-5 regenerated but 0 PASS records "
                               "(old design?)")

    # Legacy outputs present but not used by the current chain
    legacy = sorted((data["dir"] / "4.SNV").glob("*.norm.snpeff.vcf"))
    if legacy and not data["filtering"]["mutect2_pass"]:
        for path in legacy[:1]:
            n = None
            try:
                with open(path, errors="replace") as fh:
                    n = sum(1 for line in fh if not line.startswith("#"))
            except OSError:
                pass
            data["notes"].append(
                f"Legacy snpEff output present (NOT current chain): "
                f"{path.relative_to(data['dir'])} ({fmt_int(n)} records)")

    status["hla"] = ("ok", "") if data["hla"].get("n_loci") else (
        "missing", "R1_bestguess_G.txt missing")

    return status


def overall_status(status: dict) -> str:
    rank = {"ok": 0, "incomplete": 1, "stale": 2, "missing": 3}
    worst = "ok"
    for value in status.values():
        if isinstance(value, tuple):
            if rank.get(value[0], 9) > rank.get(worst, 0):
                worst = value[0]
        else:
            for v in value.values():
                if rank.get(v[0], 9) > rank.get(worst, 0):
                    worst = v[0]
    return worst


STATUS_BADGE = {"ok": "[OK]", "incomplete": "[INCOMPLETE]",
                "stale": "[STALE]", "missing": "[MISSING]"}


def badge(status: dict, key: str) -> str:
    val = status.get(key)
    code = val[0] if isinstance(val, tuple) else "ok"
    return STATUS_BADGE.get(code, "")


# ---------------------------------------------------------------------------
# Renderers
# ---------------------------------------------------------------------------

def render_patient_md(data: dict) -> str:
    p = data["patient"]
    st = data["status"]
    out: list[str] = [f"# {p}", ""]

    # Data status
    out += ["## Data status", ""]
    rows = []
    for key, label in [("qc", "1. QC (fastp)"), ("alignment", "2. Alignment (BWA-MEM2)"),
                       ("bqsr", "3. BQSR"), ("contamination", "3.2 Contamination"),
                       ("calling", "4. SNV calling"), ("filtering", "5. Filtering"),
                       ("annotation", "5. Annotation (snpEff)"), ("hla", "7. HLA typing")]:
        val = st.get(key)
        if isinstance(val, tuple):
            code, reason = val
        else:
            codes = {v[0] for v in val.values()}
            code = "ok" if codes == {"ok"} else sorted(codes)[0]
            reason = "; ".join(v[1] for v in val.values() if v[1])
        rows.append([label, STATUS_BADGE.get(code, code), reason or "-"])
    out += [md_table(["Stage", "Status", "Note"], rows), ""]
    if data["notes"]:
        out += ["> **Notes**", ">"] + [f"> - {n}" for n in data["notes"]] + [""]

    # 1. QC
    out += [f"## 1. QC (fastp) {badge(st, 'qc')}", ""]
    qc_rows = []
    for side in ("tumor", "blood"):
        s = data["sides"][side]
        q = s["qc"]
        qc_rows.append([
            side,
            fmt_int(q.get("raw_reads")), fmt_int(q.get("clean_reads")),
            fmt_pct(q.get("pass_rate")),
            fmt_pct(q.get("q30_before")), fmt_pct(q.get("q30_after")),
            fmt_pct(q.get("gc_after")),
            fmt_pct(q.get("dup_rate")),
            fmt_int(q.get("adapter_trimmed_reads")),
            fmt_int(q.get("insert_peak")),
        ])
    out += [md_table(["Side", "Raw reads", "Clean reads", "Pass rate", "Q30 before",
                      "Q30 after", "GC after", "Dup rate", "Adapter-trimmed reads",
                      "Insert peak"], qc_rows), ""]
    for side in ("tumor", "blood"):
        q = data["sides"][side]["qc"]
        if q.get("passed_filter_reads") is not None:
            out.append(f"- **{side}** filtered reads: "
                       f"passed {fmt_int(q.get('passed_filter_reads'))}, "
                       f"low-quality {fmt_int(q.get('low_quality_reads'))}, "
                       f"too-short {fmt_int(q.get('too_short_reads'))}, "
                       f"adapter-dimer {fmt_int(q.get('adapter_dimer_reads'))}")
    out.append("")

    # 2. Alignment
    out += [f"## 2. Alignment and BQSR {badge(st, 'alignment')}", ""]
    al_rows = []
    for side in ("tumor", "blood"):
        s = data["sides"][side]
        for label, fs in [("BWA raw", s["raw_flagstat"]),
                          ("markdup", s["markdup_flagstat"]),
                          ("BQSR", s["recal_flagstat"])]:
            total = fs.get("in total", {}).get("count")
            mapped = fs.get("mapped", {})
            prop = fs.get("properly paired", {})
            dup = fs.get("duplicates", {})
            al_rows.append([
                side, label,
                fmt_int(total), fmt_int(mapped.get("count")),
                f"{mapped['pct']:.2f}%" if mapped.get("pct") is not None else NULL_MARK,
                f"{prop['pct']:.2f}%" if prop.get("pct") is not None else NULL_MARK,
                fmt_int(dup.get("count")),
                fmt_int(fs.get("singletons", {}).get("count")),
            ])
    out += [md_table(["Side", "BAM", "Total", "Mapped", "Mapped %",
                      "Properly paired %", "Duplicates", "Singletons"], al_rows), ""]
    pic_rows = []
    for side in ("tumor", "blood"):
        s = data["sides"][side]
        pic_rows.append([side, fmt_pct(s["dup_pct"]),
                         fmt_int(safe_int(s["picard"].get("ESTIMATED_LIBRARY_SIZE")))])
    out += [md_table(["Side", "Picard duplication", "Est. library size"], pic_rows), ""]

    # 3. Contamination
    out += [f"## 3. Contamination {badge(st, 'contamination')}", ""]
    out += [f"- Tumor contamination: **{fmt_float(data['contamination']['contamination'], 6)}**",
            f"- Error rate: {fmt_float(data['contamination']['error'], 8)}", ""]

    # 4. Calling
    out += [f"## 4. SNV calling {badge(st, 'calling')}", ""]
    c = data["calling"]
    strl = c["strelka"]
    snv_tally = strl["snvs"]["filters"]
    pass_snv = snv_tally.get("PASS", 0)
    out += [f"- Mutect2 unfiltered sites: **{fmt_int(c['unfiltered'])}**"]
    out += [f"- Strelka2 SNV sites: **{fmt_int(strl['snvs']['n'])}** "
            f"(PASS {fmt_int(pass_snv)})"]
    out += [f"- Strelka2 indel sites: **{fmt_int(strl['indels']['n'])}** "
            f"(PASS {fmt_int(strl['indels']['filters'].get('PASS', 0))})"]
    out += ["", "Strelka2 SNV FILTER breakdown:", ""]
    out += [md_table(["FILTER", "Records"],
                     [[k, fmt_int(v)] for k, v in snv_tally.most_common()]), ""]

    # 5. Filtering
    filt = data["filtering"]
    out += [f"## 5. Filtering {badge(st, 'filtering')}", ""]
    meta = filt["stats"].get("meta", {})
    out += [f"- Mutect2 PASS sites: **{fmt_int(filt['mutect2_pass'])}**",
            f"- Both Mutect2 + Strelka2 PASS: **{fmt_int(filt['both_pass'])}**",
            f"- Filter threshold: {fmt_float(meta.get('threshold'), 3)}, "
            f"FDR: {fmt_float(meta.get('fdr'), 3)}, "
            f"sensitivity: {fmt_float(meta.get('sensitivity'), 3)}", ""]
    if filt["stats"]["filters"]:
        out += [md_table(["Filter", "FP", "FDR", "FN", "FNR"],
                         [[f["filter"], fmt_float(f["FP"], 2), fmt_float(f["FDR"], 2),
                           fmt_float(f["FN"], 2), fmt_float(f["FNR"], 2)]
                          for f in filt["stats"]["filters"]]), ""]

    # 6. Annotation
    out += [f"## 6. Annotation (snpEff) {badge(st, 'annotation')}", ""]
    ann_rows = []
    for chain, label in CHAIN_LABELS.items():
        info = data["chains"][chain]
        mf, pd = info["missense_frameshift"], info["protein_damaging"]
        ann_rows.append([
            label,
            fmt_int(mf["n_sites"]) if mf["readable"] else NULL_MARK,
            fmt_int(pd["n_sites"]) if pd["readable"] else NULL_MARK,
            fmt_float(median(pd["tlods"]), 2) if pd["readable"] else NULL_MARK,
            fmt_int(pd["impacts"].get("HIGH", 0)) if pd["readable"] else NULL_MARK,
            fmt_int(pd["impacts"].get("MODERATE", 0)) if pd["readable"] else NULL_MARK,
        ])
    out += [md_table(["Chain", "Missense/frameshift sites", "Protein-damaging sites",
                      "Median TLOD", "HIGH (transcripts)",
                      "MODERATE (transcripts)"], ann_rows), ""]
    for chain, label in CHAIN_LABELS.items():
        pd = data["chains"][chain]["protein_damaging"]
        if pd["readable"] and pd["effect_any"]:
            top = ", ".join(f"{k} ({v})" for k, v in pd["effect_any"].most_common(8))
            out += [f"- **{label}** top effects: {top}"]
    out.append("")

    # 7. Cancer genes
    out += [f"## 7. Cancer-associated genes", ""]
    total_genes = len(data["genes"])
    flagged = [g for g in data["genes"] if g["cancer"]]
    by_tier = Counter(g["cancer"]["tier"] for g in flagged)
    out += [f"Gene source: `{data['gene_source']}` — {fmt_int(total_genes)} genes hit, "
            f"**{fmt_int(len(flagged))}** present in at least one cancer gene list "
            f"(tier 3: {by_tier.get(3, 0)}, tier 2: {by_tier.get(2, 0)}, "
            f"tier 1: {by_tier.get(1, 0)}).", ""]
    if data["genes"]:
        rows = []
        for g in data["genes"]:
            rows.append([
                g["gene"], fmt_int(g["n_sites"]), fmt_int(g["high_sites"]),
                g["worst_impact"], g["worst_consequence"], fmt_float(g["tlod_max"], 2),
                cancer_cell(g["cancer"]),
            ])
        out += [md_table(["Gene", "Sites", "HIGH sites", "Worst impact",
                          "Worst consequence", "Max TLOD", "Cancer gene"], rows), ""]
    else:
        out += ["_No gene-level data available._", ""]

    # 8. HLA
    hla = data["hla"]
    out += [f"## 8. HLA typing {badge(st, 'hla')}", ""]
    if hla.get("n_loci"):
        rows = []
        for locus in hla["loci"]:
            alleles = " / ".join(";".join(a) if a else NULL_MARK for a in locus["alleles"])
            cov = [c for c in locus["coverage"] if c is not None]
            q1 = [q for q in locus["q1"] if q is not None]
            rows.append([
                locus["locus"], alleles,
                fmt_float(q1[0], 3) if q1 else NULL_MARK,
                fmt_float(cov[0], 1) if cov else NULL_MARK,
                f"{locus['n_perfect']}/{len(locus['perfect'])}",
                "yes" if locus["ambiguous"] else "",
            ])
        out += [md_table(["Locus", "Alleles (hap1 / hap2)", "Q1", "Avg coverage",
                          "perfectG", "Ambiguous"], rows), ""]
        out += [f"- Loci typed: **{hla['n_loci']}**, perfectG loci: "
                f"**{hla['n_perfect_g']}**, ambiguous loci: {hla['n_ambiguous']}", ""]
    else:
        out += ["_No HLA results found._", ""]

    return "\n".join(out)


def gene_table_tsv(data: dict) -> str:
    lines = ["gene\tn_sites\thigh_sites\tworst_impact\tworst_consequence\ttlod_max\tcancer_sources"]
    for g in data["genes"]:
        srcs = ",".join(g["cancer"]["sources"]) if g["cancer"] else ""
        lines.append("\t".join([
            g["gene"], str(g["n_sites"]), str(g["high_sites"]),
            g["worst_impact"], g["worst_consequence"],
            "" if g["tlod_max"] is None else f"{g['tlod_max']:.4g}", srcs,
        ]))
    return "\n".join(lines) + "\n"


COHORT_COLUMNS = [
    "patient", "status", "notes",
    "tumor_raw_reads", "tumor_clean_reads", "tumor_pass_rate", "tumor_q30_after",
    "tumor_gc_after", "blood_raw_reads", "blood_clean_reads", "blood_pass_rate",
    "tumor_mapped_pct", "tumor_dup_pct", "tumor_properly_paired_pct",
    "blood_mapped_pct", "blood_dup_pct",
    "contamination", "unfiltered_sites", "strelka_snv_sites", "strelka_pass",
    "mutect2_pass", "both_pass", "fdr", "sensitivity",
    "missense_frameshift_mutect2", "protein_damaging_mutect2",
    "missense_frameshift_both", "protein_damaging_both",
    "genes_hit_mutect2", "cancer_genes_mutect2", "cancer_tier3",
    "cancer_tier2", "cancer_tier1", "median_tlod_mutect2",
    "hla_loci", "hla_perfect_g_loci",
]


def cohort_row(data: dict) -> dict:
    def side(s):
        return data["sides"][s]

    def q(s, k):
        return side(s)["qc"].get(k)

    def fs(s, key, sub):
        return side(s)["markdup_flagstat"].get(key, {}).get(sub)

    filt = data["filtering"]
    chains = data["chains"]
    flagged = [g for g in data["genes"] if g["cancer"]]
    tiers = Counter(g["cancer"]["tier"] for g in flagged)
    pd = chains["mutect2"]["protein_damaging"]

    return {
        "patient": data["patient"],
        "status": overall_status(data["status"]),
        "notes": " | ".join(data["notes"]),
        "tumor_raw_reads": q("tumor", "raw_reads"),
        "tumor_clean_reads": q("tumor", "clean_reads"),
        "tumor_pass_rate": round(side("tumor")["qc"].get("pass_rate"), 6)
        if side("tumor")["qc"].get("pass_rate") is not None else None,
        "tumor_q30_after": q("tumor", "q30_after"),
        "tumor_gc_after": q("tumor", "gc_after"),
        "blood_raw_reads": q("blood", "raw_reads"),
        "blood_clean_reads": q("blood", "clean_reads"),
        "blood_pass_rate": round(side("blood")["qc"].get("pass_rate"), 6)
        if side("blood")["qc"].get("pass_rate") is not None else None,
        "tumor_mapped_pct": fs("tumor", "mapped", "pct"),
        "tumor_dup_pct": side("tumor")["dup_pct"],
        "tumor_properly_paired_pct": fs("tumor", "properly paired", "pct"),
        "blood_mapped_pct": fs("blood", "mapped", "pct"),
        "blood_dup_pct": side("blood")["dup_pct"],
        "contamination": data["contamination"]["contamination"],
        "unfiltered_sites": data["calling"]["unfiltered"],
        "strelka_snv_sites": data["calling"]["strelka"]["snvs"]["n"],
        "strelka_pass": data["calling"]["strelka"]["snvs"]["filters"].get("PASS"),
        "mutect2_pass": filt["mutect2_pass"],
        "both_pass": filt["both_pass"],
        "fdr": filt["stats"]["meta"].get("fdr"),
        "sensitivity": filt["stats"]["meta"].get("sensitivity"),
        "missense_frameshift_mutect2": chains["mutect2"]["missense_frameshift"]["n_sites"],
        "protein_damaging_mutect2": pd["n_sites"] if pd["readable"] else None,
        "missense_frameshift_both": chains["both"]["missense_frameshift"]["n_sites"],
        "protein_damaging_both": chains["both"]["protein_damaging"]["n_sites"]
        if chains["both"]["protein_damaging"]["readable"] else None,
        "genes_hit_mutect2": len(data["genes"]) or None,
        "cancer_genes_mutect2": len(flagged) if data["genes"] else None,
        "cancer_tier3": tiers.get(3, 0) if data["genes"] else None,
        "cancer_tier2": tiers.get(2, 0) if data["genes"] else None,
        "cancer_tier1": tiers.get(1, 0) if data["genes"] else None,
        "median_tlod_mutect2": median(pd["tlods"]) if pd["readable"] else None,
        "hla_loci": data["hla"].get("n_loci"),
        "hla_perfect_g_loci": data["hla"].get("n_perfect_g"),
    }


def render_summary_md(rows: list[dict], gene_stats: dict, patients: list[str]) -> str:
    out = ["# OE_WGS cohort summary", "",
           f"Patients: **{len(patients)}** — {', '.join(patients)}", ""]

    out += ["## Cancer gene lists", ""]
    src = gene_stats.get("sources", {})
    out += [f"- Union across the three sources: **{fmt_int(gene_stats.get('union'))}** genes",
            f"- Present in all three: **{fmt_int(gene_stats.get('all_three'))}**",
            f"- OncoKB: {fmt_int(src.get('OncoKB'))}, IntOGen: {fmt_int(src.get('IntOGen'))}, "
            f"COSMIC: {fmt_int(src.get('COSMIC'))}", ""]
    for f in gene_stats.get("files", []):
        out.append(f"- `{f}`")
    out += [""]

    def col(key, label, fmt=fmt_int):
        return [label, [fmt(r.get(key)) for r in rows]]

    out += ["## QC and alignment", ""]
    out += [md_table(
        ["patient", "status", "raw reads (T)", "clean reads (T)", "pass % (T)",
         "Q30 after (T)", "mapped % (T)", "dup % (T)", "proper pair % (T)",
         "raw reads (B)", "mapped % (B)"],
        [[r["patient"], r["status"], fmt_int(r["tumor_raw_reads"]),
          fmt_int(r["tumor_clean_reads"]), fmt_pct(r["tumor_pass_rate"]),
          fmt_pct(r["tumor_q30_after"]),
          f"{r['tumor_mapped_pct']:.2f}%" if r["tumor_mapped_pct"] is not None else NULL_MARK,
          fmt_pct(r["tumor_dup_pct"]),
          f"{r['tumor_properly_paired_pct']:.2f}%"
          if r["tumor_properly_paired_pct"] is not None else NULL_MARK,
          fmt_int(r["blood_raw_reads"]),
          f"{r['blood_mapped_pct']:.2f}%" if r["blood_mapped_pct"] is not None else NULL_MARK]
         for r in rows]), ""]

    out += ["## Variant calling and filtering", ""]
    out += [md_table(
        ["patient", "unfiltered", "Strelka SNV", "Strelka PASS", "Mutect2 PASS",
         "both PASS", "FDR", "sensitivity", "contamination"],
        [[r["patient"], fmt_int(r["unfiltered_sites"]), fmt_int(r["strelka_snv_sites"]),
          fmt_int(r["strelka_pass"]), fmt_int(r["mutect2_pass"]), fmt_int(r["both_pass"]),
          fmt_float(r["fdr"], 3), fmt_float(r["sensitivity"], 3),
          fmt_float(r["contamination"], 6)] for r in rows]), ""]

    out += ["## Annotation and cancer genes", ""]
    out += [md_table(
        ["patient", "missense/fs (M)", "protein-damaging (M)", "missense/fs (both)",
         "protein-damaging (both)", "genes hit", "cancer genes", "tier3", "tier2",
         "tier1", "median TLOD"],
        [[r["patient"], fmt_int(r["missense_frameshift_mutect2"]),
          fmt_int(r["protein_damaging_mutect2"]),
          fmt_int(r["missense_frameshift_both"]), fmt_int(r["protein_damaging_both"]),
          fmt_int(r["genes_hit_mutect2"]), fmt_int(r["cancer_genes_mutect2"]),
          fmt_int(r["cancer_tier3"]), fmt_int(r["cancer_tier2"]),
          fmt_int(r["cancer_tier1"]), fmt_float(r["median_tlod_mutect2"], 2)]
         for r in rows]), ""]

    out += ["## HLA typing", ""]
    out += [md_table(["patient", "loci typed", "perfectG loci"],
                     [[r["patient"], fmt_int(r["hla_loci"]),
                       fmt_int(r["hla_perfect_g_loci"])] for r in rows]), ""]

    flagged = [r for r in rows if r["status"] != "ok"]
    if flagged:
        out += ["## Patients needing attention", ""]
        out += [md_table(["patient", "status", "notes"],
                         [[r["patient"], r["status"], r["notes"] or "-"]
                          for r in flagged]), ""]
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def run_verify(data: dict, method: str, bcftools: str | None, strict: bool) -> int:
    failures = 0

    def check(name, ok, detail=""):
        nonlocal failures
        if ok:
            print(f"VERIFY ok   {name} {detail}")
        else:
            failures += 1
            print(f"VERIFY FAIL {name} {detail}")

    p = data["patient"]
    print(f"--- verifying {p} ---")

    # 1. counting methods agree
    if method == "auto" and bcftools:
        for label, path in [
            ("unfiltered", data["calling"]["unfiltered_path"]),
            ("mutect2_pass", data["filtering"]["mutect2_pass_path"]),
            ("both_pass", data["filtering"]["both_pass_path"]),
            ("strelka_snv", data["calling"]["strelka"]["snvs"]["path"]),
        ]:
            if path is None or not path.exists():
                continue
            a = count_vcf_records(path, "bcftools", bcftools)
            b = count_vcf_records(path, "gzip", None)
            check(f"count[{label}]", a == b, f"bcftools={a} gzip={b}")

    # 2. fastp internal consistency
    for side in ("tumor", "blood"):
        q = data["sides"][side]["qc"]
        if q.get("passed_filter_reads") is not None and q.get("clean_reads") is not None:
            check(f"fastp[{side}] passed==clean",
                  q["passed_filter_reads"] == q["clean_reads"],
                  f"{q['passed_filter_reads']} vs {q['clean_reads']}")
            check(f"fastp[{side}] after<before",
                  q["clean_reads"] < q["raw_reads"],
                  f"{q['clean_reads']} < {q['raw_reads']}")

    # 3. duplication cross-check
    for side in ("tumor", "blood"):
        s = data["sides"][side]
        picard = safe_float(s["picard"].get("PERCENT_DUPLICATION"))
        fs = s["markdup_flagstat"]
        total = fs.get("in total", {}).get("count")
        dup = fs.get("duplicates", {}).get("count")
        if picard is not None and total:
            derived = dup / total
            check(f"dup[{side}] picard~=flagstat", abs(picard - derived) < 0.005,
                  f"picard={picard:.5f} flagstat={derived:.5f}")

    # 4. gene table integrity
    for chain in ("mutect2", "both"):
        pd = data["chains"][chain]["protein_damaging"]
        if pd["readable"]:
            check(f"genes[{chain}] header", pd["header_ok"])
            sites = {s for e in pd["genes"].values() for s in e["sites"]}
            check(f"genes[{chain}] sites<=n_sites", len(sites) <= pd["n_sites"],
                  f"{len(sites)} gene-linked of {pd['n_sites']}")
            if pd["n_sites"]:
                check(f"genes[{chain}] n_genes", len(pd["genes"]) > 0,
                      f"{len(pd['genes'])} genes")
            else:
                # Header-only TSV is a valid "nothing passed" result.
                check(f"genes[{chain}] empty-consistent", not pd["genes"],
                      "0 sites -> 0 genes")

    # 5. HLA shape
    hla = data["hla"]
    if hla.get("n_loci"):
        check("hla loci==13", hla["n_loci"] == 13, f"{hla['n_loci']} loci")
        bad = [l for l in hla["loci"] if len(l["alleles"]) > 2]
        check("hla <=2 alleles/locus", not bad, f"{len(bad)} loci with >2")

    print(f"--- {failures} failure(s) ---")
    return failures


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def discover_patients(base: Path) -> list[str]:
    out = []
    for entry in sorted(base.iterdir()):
        if not entry.is_dir() or entry.name in SKIP_DIRS or entry.name.startswith("."):
            continue
        if any((entry / s).is_dir() for s in ("WGS", "WGS_blood", "5.FilterSNV", "7.HLA_type")):
            out.append(entry.name)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", default=".", type=Path)
    ap.add_argument("--patients", nargs="*", default=None)
    ap.add_argument("--outdir", default="reports", type=Path)
    ap.add_argument("--genes-dir", default="CancerGens", type=Path)
    ap.add_argument("--count-method", choices=["auto", "bcftools", "gzip"], default="auto")
    ap.add_argument("--bcftools", default=None)
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args(argv)

    base = args.base_dir.resolve()
    genes_dir = (base / args.genes_dir) if not args.genes_dir.is_absolute() else args.genes_dir
    outdir = (base / args.outdir) if not args.outdir.is_absolute() else args.outdir
    tables_dir = outdir / "tables"

    bcftools = args.bcftools or shutil.which("bcftools")
    method = args.count_method
    if method in ("auto", "bcftools") and not bcftools:
        if method == "bcftools":
            print("warning: bcftools not found, falling back to gzip", file=sys.stderr)
        method = "gzip"

    patients = args.patients or discover_patients(base)
    if not patients:
        print("error: no patient directories found", file=sys.stderr)
        return 1

    genes, gene_stats = load_cancer_genes(genes_dir)
    if not genes:
        print(f"warning: no cancer gene lists loaded from {genes_dir}", file=sys.stderr)
    else:
        print(f"loaded {gene_stats['union']} cancer genes "
              f"({gene_stats['all_three']} in all three sources) from {genes_dir}",
              file=sys.stderr)

    outdir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    all_data = []
    for patient in patients:
        pdir = base / patient
        if not pdir.is_dir():
            print(f"warning: {pdir} not found, skipping", file=sys.stderr)
            continue
        print(f"extracting {patient} ...", file=sys.stderr)
        data = extract_patient(patient, base, genes, method, bcftools)
        all_data.append(data)

        (outdir / f"{patient}.report.md").write_text(render_patient_md(data))
        (tables_dir / f"{patient}.genes.tsv").write_text(gene_table_tsv(data))

    rows = [cohort_row(d) for d in all_data]
    with open(outdir / "cohort_summary.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COHORT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k))
                             for k in COHORT_COLUMNS})

    (outdir / "SUMMARY.md").write_text(
        render_summary_md(rows, gene_stats, [d["patient"] for d in all_data]))

    print(f"\nwrote {len(all_data)} patient reports + cohort_summary.tsv + SUMMARY.md "
          f"to {outdir}", file=sys.stderr)

    if args.verify:
        total = 0
        for data in all_data:
            total += run_verify(data, method, bcftools, args.strict)
        if args.strict and total:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
