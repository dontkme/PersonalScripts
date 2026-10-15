#!/bin/bash
# Link outputs from one pipeline stage into the next stage's directory.
#
# Usage:
#   bash link_stage_outputs.sh <sample|all> <stage>
#
#   stage 1:   clean.fq.gz        from 1.QC/          -> 2.BWA2_res/
#   stage 2:   markdup.bam/.bai   from 2.BWA2_res/    -> 3.BQSR/
#   stage 3:   recal.bam/.bai     from 3.BQSR/        -> 3.1Pileup/, 4.SNV/, 4_1.St_SNV/, 6.CNV/, 7.HLA_type/ (blood)
#   stage 3.1: pileups.table      from 3.1Pileup*/    -> 3.2CalCont/
#   stage 3.2: contamination/segs from 3.2CalCont/    -> 5.FilterSNV/
#   stage 4:   vcf.gz/.stats/.tbi from 4.SNV/ + Strelka2 VCFs from 4_1.St_SNV/ -> 5.FilterSNV/
#   stage 4.f1r2: read-orientation from 4.SNV/        -> 5.FilterSNV/
#   stage 6.1:  FACETS VCF (6.CNV/) + SSM VCF (5.FilterSNV/) -> 6.1phyloWGS/
#   stage 8:    both-PASS VCF (5.FilterSNV/) + FACETS VCF (6.CNV/) + HLA (7.HLA_type/) + phyloWGS top trees (6.1phyloWGS/) -> 8.VacQ/

set -euo pipefail

BASE_DIR="/sibcb1/jizhelab1/hukaining/202608/OE_WGS"
#ALL_SAMPLES=(B1C391 BDZ717 BEC478 BEF141 BEG091 BEG323 BEH047 E89007)
ALL_SAMPLES=(BE0728 E88354 E88052 B1C391 BDZ717 BEC478 BEF141 BEG091 BEG323 BEH047 E89007)

die() { echo "ERROR: $*" >&2; exit 1; }

if [[ $# -lt 2 ]]; then
    echo "Usage: bash link_stage_outputs.sh <sample|all> <stage>"
    echo ""
    echo "  stage 1:   clean.fq.gz        from 1.QC/          -> 2.BWA2_res/"
    echo "  stage 2:   markdup.bam/.bai   from 2.BWA2_res/    -> 3.BQSR/"
    echo "  stage 3:   recal.bam/.bai     from 3.BQSR/        -> 3.1Pileup/, 4.SNV/, 4_1.St_SNV/, 6.CNV/, 7.HLA_type/"
    echo "  stage 3.1: pileups.table      from 3.1Pileup*/    -> 3.2CalCont/"
    echo "  stage 3.2: contamination/segs from 3.2CalCont/    -> 5.FilterSNV/"
    echo "  stage 4:   vcf.gz/.stats/.tbi from 4.SNV/ + Strelka2 VCFs from 4_1.St_SNV/ -> 5.FilterSNV/"
    echo "  stage 4.f1r2: read-orientation  from 4.SNV/         -> 5.FilterSNV/"
    echo "  stage 6.1:  FACETS VCF + SSM VCF                   -> 6.1phyloWGS/"
    echo "  stage 8:    both-PASS VCF + FACETS VCF + HLA       -> 8.VacQ/"
    exit 1
fi

SAMPLE_INPUT="$1"
STAGE="$2"

if [[ "$SAMPLE_INPUT" == "all" ]]; then
    SAMPLES=("${ALL_SAMPLES[@]}")
else
    SAMPLES=("$SAMPLE_INPUT")
fi

for sample in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${BASE_DIR}/${sample}"

    if [[ ! -d "$SAMPLE_DIR" ]]; then
        echo "WARNING: ${SAMPLE_DIR} does not exist, skipping"
        continue
    fi

    echo "--- Linking stage ${STAGE} outputs for ${sample} ---"

    # ==================================================================
    # Stage 1: fastp outputs -> 2.BWA2_res
    # ==================================================================
    if [[ "$STAGE" == "1" ]]; then
        for sub in WGS WGS_blood; do
            SRC="${SAMPLE_DIR}/${sub}/1.QC"
            DST="${SAMPLE_DIR}/${sub}/2.BWA2_res"
            if [[ -d "$SRC" ]] && [[ -d "$DST" ]]; then
                for f in "$SRC"/*.clean.fq.gz; do
                    if [[ -f "$f" ]]; then
                        base=$(basename "$f")
                        ln -sf "../1.QC/${base}" "${DST}/${base}"
                        echo "  ${sub}: ${base} -> 2.BWA2_res/"
                    fi
                done
            fi
        done

    # ==================================================================
    # Stage 2: BWA2/markdup outputs -> 3.BQSR
    # ==================================================================
    elif [[ "$STAGE" == "2" ]]; then
        for sub in WGS WGS_blood; do
            SRC="${SAMPLE_DIR}/${sub}/2.BWA2_res"
            DST="${SAMPLE_DIR}/${sub}/3.BQSR"
            if [[ -d "$SRC" ]] && [[ -d "$DST" ]]; then
                for f in "$SRC"/*.markdup.bam "$SRC"/*.markdup.bai; do
                    if [[ -f "$f" ]]; then
                        base=$(basename "$f")
                        ln -sf "../2.BWA2_res/${base}" "${DST}/${base}"
                        echo "  ${sub}: ${base} -> 3.BQSR/"
                    fi
                done
            fi
        done

    # ==================================================================
    # Stage 3: BQSR outputs -> 3.1Pileup, 4.SNV, 4_1.St_SNV
    # ==================================================================
    elif [[ "$STAGE" == "3" ]]; then
        # WGS tumor side: 3.BQSR -> WGS/3.1Pileup_Conta
        SRC="${SAMPLE_DIR}/WGS/3.BQSR"
        for dst_rel in "WGS/3.1Pileup_Conta"; do
            DST="${SAMPLE_DIR}/${dst_rel}"
            if [[ -d "$SRC" ]] && [[ -d "$DST" ]]; then
                for f in "$SRC"/*.markdup.recal.bam "$SRC"/*.markdup.recal.bai "$SRC"/*.markdup.recal.table "$SRC"/*.markdup.recal.bam.flagstat; do
                    if [[ -f "$f" ]]; then
                        base=$(basename "$f")
                        # WGS/3.1Pileup_Conta links to ../3.BQSR/ (it's in WGS/)
                        ln -sf "../3.BQSR/${base}" "${DST}/${base}"
                        echo "  ${dst_rel}: ${base}"
                    fi
                done
            fi
        done

        # WGS_blood side: 3.BQSR -> WGS_blood/3.1Pileup
        SRC_B="${SAMPLE_DIR}/WGS_blood/3.BQSR"
        DST_B="${SAMPLE_DIR}/WGS_blood/3.1Pileup"
        if [[ -d "$SRC_B" ]] && [[ -d "$DST_B" ]]; then
            for f in "$SRC_B"/*.markdup.recal.bam "$SRC_B"/*.markdup.recal.bai "$SRC_B"/*.markdup.recal.table "$SRC_B"/*.markdup.recal.bam.flagstat; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../3.BQSR/${base}" "${DST_B}/${base}"
                    echo "  WGS_blood/3.1Pileup: ${base}"
                fi
            done
        fi

        # Top-level stages: link from both WGS/3.BQSR and WGS_blood/3.BQSR
        # (6.CNV needs both BAMs for FACETS allele-specific copy number)
        for top_stage in "4.SNV" "4_1.St_SNV" "6.CNV"; do
            DST_TOP="${SAMPLE_DIR}/${top_stage}"
            if [[ ! -d "$DST_TOP" ]]; then continue; fi
            # Tumor recal files
            if [[ -d "$SRC" ]]; then
                for f in "$SRC"/*.markdup.recal.bam "$SRC"/*.markdup.recal.bai "$SRC"/*.markdup.recal.table "$SRC"/*.markdup.recal.bam.flagstat; do
                    if [[ -f "$f" ]]; then
                        base=$(basename "$f")
                        ln -sf "../WGS/3.BQSR/${base}" "${DST_TOP}/${base}"
                        echo "  ${top_stage}: ${base} (tumor)"
                    fi
                done
            fi
            # Blood recal files
            if [[ -d "$SRC_B" ]]; then
                for f in "$SRC_B"/*.markdup.recal.bam "$SRC_B"/*.markdup.recal.bai "$SRC_B"/*.markdup.recal.table "$SRC_B"/*.markdup.recal.bam.flagstat; do
                    if [[ -f "$f" ]]; then
                        base=$(basename "$f")
                        ln -sf "../WGS_blood/3.BQSR/${base}" "${DST_TOP}/${base}"
                        echo "  ${top_stage}: ${base} (blood)"
                    fi
                done
            fi
        done

        # 7.HLA_type: link the blood recal BAM only (HLA typing is germline)
        DST_HLA="${SAMPLE_DIR}/7.HLA_type"
        if [[ -d "$SRC_B" ]] && [[ -d "$DST_HLA" ]]; then
            for f in "$SRC_B"/*.markdup.recal.bam "$SRC_B"/*.markdup.recal.bai "$SRC_B"/*.markdup.recal.bam.flagstat; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../WGS_blood/3.BQSR/${base}" "${DST_HLA}/${base}"
                    echo "  7.HLA_type: ${base} (blood)"
                fi
            done
        fi

    # ==================================================================
    # Stage 3.1: Pileup outputs -> 3.2CalCont
    # ==================================================================
    elif [[ "$STAGE" == "3.1" ]]; then
        DST="${SAMPLE_DIR}/3.2CalCont"
        if [[ ! -d "$DST" ]]; then
            die "3.2CalCont directory not found for ${sample}"
        fi

        # Tumor pileup
        SRC_T="${SAMPLE_DIR}/WGS/3.1Pileup_Conta"
        if [[ -d "$SRC_T" ]]; then
            for f in "$SRC_T"/*.pileups.table; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../WGS/3.1Pileup_Conta/${base}" "${DST}/${base}"
                    echo "  Tumor: ${base} -> 3.2CalCont/"
                fi
            done
        fi

        # Blood pileup
        SRC_B="${SAMPLE_DIR}/WGS_blood/3.1Pileup"
        if [[ -d "$SRC_B" ]]; then
            for f in "$SRC_B"/*.pileups.table; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../WGS_blood/3.1Pileup/${base}" "${DST}/${base}"
                    echo "  Blood: ${base} -> 3.2CalCont/"
                fi
            done
        fi

    # ==================================================================
    # Stage 3.2: CalCont outputs -> 5.FilterSNV
    # ==================================================================
    elif [[ "$STAGE" == "3.2" ]]; then
        SRC="${SAMPLE_DIR}/3.2CalCont"
        DST="${SAMPLE_DIR}/5.FilterSNV"
        if [[ -d "$SRC" ]] && [[ -d "$DST" ]]; then
            for f in "$SRC"/*.contamination.table "$SRC"/*.segments.table; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../3.2CalCont/${base}" "${DST}/${base}"
                    echo "  ${base} -> 5.FilterSNV/"
                fi
            done
        fi

    # ==================================================================
    # Stage 4: Mutect2 outputs -> 5.FilterSNV
    # ==================================================================
    elif [[ "$STAGE" == "4" ]]; then
        SRC="${SAMPLE_DIR}/4.SNV"
        DST="${SAMPLE_DIR}/5.FilterSNV"
        if [[ -d "$SRC" ]] && [[ -d "$DST" ]]; then
            for f in "$SRC"/*.vcf.gz "$SRC"/*.vcf.gz.stats "$SRC"/*.vcf.gz.tbi; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../4.SNV/${base}" "${DST}/${base}"
                    echo "  ${base} -> 5.FilterSNV/"
                fi
            done
        fi

        # Strelka2 somatic VCFs (nested under the run dir) -> 5.FilterSNV,
        # renamed to <sample>.strelka.somatic.{snvs,indels}.vcf.gz
        ST_SRC="${SAMPLE_DIR}/4_1.St_SNV/${sample}.st_somatic_BQSR/results/variants"
        if [[ -d "$ST_SRC" ]] && [[ -d "$DST" ]]; then
            for kind in snvs indels; do
                for ext in vcf.gz vcf.gz.tbi; do
                    f="${ST_SRC}/somatic.${kind}.${ext}"
                    if [[ -f "$f" ]]; then
                        ln -sf "../4_1.St_SNV/${sample}.st_somatic_BQSR/results/variants/somatic.${kind}.${ext}" \
                               "${DST}/${sample}.strelka.somatic.${kind}.${ext}"
                        echo "  ${sample}.strelka.somatic.${kind}.${ext} -> 5.FilterSNV/ (Strelka2)"
                    fi
                done
            done
        fi

    # ==================================================================
    # Stage 4.f1r2: LearnReadOrientationModel output -> 5.FilterSNV
    # ==================================================================
    elif [[ "$STAGE" == "4.f1r2" ]]; then
        SRC="${SAMPLE_DIR}/4.SNV"
        DST="${SAMPLE_DIR}/5.FilterSNV"
        if [[ -d "$SRC" ]] && [[ -d "$DST" ]]; then
            for f in "$SRC"/*.read-orientation-model.tar.gz; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../4.SNV/${base}" "${DST}/${base}"
                    echo "  ${base} -> 5.FilterSNV/"
                fi
            done
        fi

    # ==================================================================
    # Stage 6.1: phyloWGS inputs (FACETS VCF + SSM VCF) -> 6.1phyloWGS
    # Run only after both 5.FilterSNV and 6.CNV have completed.
    # ==================================================================
    elif [[ "$STAGE" == "6.1" ]]; then
        DST="${SAMPLE_DIR}/6.1phyloWGS"
        if [[ ! -d "$DST" ]]; then
            die "6.1phyloWGS directory not found for ${sample}"
        fi

        # FACETS segments from step 6
        SRC_CNV="${SAMPLE_DIR}/6.CNV"
        if [[ -d "$SRC_CNV" ]]; then
            for f in "$SRC_CNV/${sample}.facet_out.vcf.gz" "$SRC_CNV/${sample}.facet_out.vcf.gz.tbi"; do
                if [[ -f "$f" ]]; then
                    base=$(basename "$f")
                    ln -sf "../6.CNV/${base}" "${DST}/${base}"
                    echo "  ${base} -> 6.1phyloWGS/ (FACETS)"
                fi
            done
        fi

        # Somatic SSM VCF from step 5
        SRC_SSM="${SAMPLE_DIR}/5.FilterSNV"
        SSM="${sample}.mutect2.filtered.PASS.norm.snpeff.missense_frameshift.vcf"
        if [[ -f "${SRC_SSM}/${SSM}" ]]; then
            ln -sf "../5.FilterSNV/${SSM}" "${DST}/${SSM}"
            echo "  ${SSM} -> 6.1phyloWGS/ (SSMs)"
        fi

    # ==================================================================
    # Stage 8: VacQ inputs (both-PASS VCF + FACETS VCF + HLA + phyloWGS) -> 8.VacQ
    # Run only after 5.FilterSNV, 6.CNV, 6.1phyloWGS and 7.HLA_type have
    # completed -- the quality scoring maps variants onto the top trees.
    # ==================================================================
    elif [[ "$STAGE" == "8" ]]; then
        DST="${SAMPLE_DIR}/8.VacQ"
        if [[ ! -d "$DST" ]]; then
            die "8.VacQ directory not found for ${sample}"
        fi

        # Both-callers-PASS missense/frameshift VCF -> preferred VEP input.
        # The Mutect2-only chain is linked as well: the batcher falls back to
        # it when the both-PASS chain filters down to zero variants.
        for VCF in "${sample}.mutect2_strelka.both.PASS.snpeff.missense_frameshift.vcf" \
                   "${sample}.mutect2.filtered.PASS.norm.snpeff.missense_frameshift.vcf"; do
            if [[ -f "${SAMPLE_DIR}/5.FilterSNV/${VCF}" ]]; then
                ln -sf "../5.FilterSNV/${VCF}" "${DST}/${VCF}"
                echo "  ${VCF} -> 8.VacQ/ (VEP input)"
            fi
        done

        # FACETS segments -> purity for --tumor-purity
        for ext in vcf.gz vcf.gz.tbi; do
            f="${SAMPLE_DIR}/6.CNV/${sample}.facet_out.${ext}"
            if [[ -f "$f" ]]; then
                ln -sf "../6.CNV/${sample}.facet_out.${ext}" "${DST}/${sample}.facet_out.${ext}"
                echo "  ${sample}.facet_out.${ext} -> 8.VacQ/ (purity)"
            fi
        done

        # phyloWGS top trees -> CCF/clonality mapping for the quality scoring
        for f in "${sample}.ssm_data.txt" "${sample}.top5.summ.json.gz" \
                 "${sample}.top5.mutass.zip" "${sample}.top5_tree_ids.tsv"; do
            if [[ -f "${SAMPLE_DIR}/6.1phyloWGS/${f}" ]]; then
                ln -sf "../6.1phyloWGS/${f}" "${DST}/${f}"
                echo "  ${f} -> 8.VacQ/ (phyloWGS top trees)"
            fi
        done

        # HLA-LA best-guess genotypes -> pVACseq allele list.
        # Layout differs per patient (BE0728 nests under <ID>_Blood/<ID>/hla).
        HLA_SRC=$(find "${SAMPLE_DIR}/7.HLA_type" -path '*/hla/R1_bestguess_G.txt' 2>/dev/null | head -1)
        if [[ -n "$HLA_SRC" ]]; then
            # 8.VacQ -> .. is the sample dir, so link relative to that
            rel="${HLA_SRC#"${SAMPLE_DIR}/"}"
            ln -sf "../${rel}" "${DST}/R1_bestguess_G.txt"
            echo "  R1_bestguess_G.txt -> 8.VacQ/ (HLA, from ${rel})"
        else
            echo "  WARNING: no hla/R1_bestguess_G.txt found for ${sample}"
        fi

    else
        die "Invalid stage: ${STAGE}. Valid: 1, 2, 3, 3.1, 3.2, 4, 4.f1r2, 6.1, 8"
    fi

    echo ""
done

echo "=== Done ==="
