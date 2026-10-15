#!/bin/bash
# One-command orchestrator for the OE_WGS somatic variant calling pipeline.
#
#   bash run_pipeline.sh <SAMPLE_ID> [options]
#
# Scaffolds the patient tree (if needed) and drives every stage to completion:
#   QC -> BWA-MEM2 -> BQSR -> pileups -> CalculateContamination
#      -> Mutect2 | Strelka2 | FACETS | HLA-LA
#      -> FilterMutectCalls + snpEff -> phyloWGS -> VEP + pVACseq + quality
#
# Each stage is submitted with `sbatch --wait --parsable`, so the driver blocks
# until it finishes and learns its true exit code (sacct is unavailable on this
# cluster). Inter-stage symlinks are created with link_stage_outputs.sh between
# stages. Run it inside tmux -- a full sample takes ~4.5 days of cluster time.
#
# Options:
#   --from <stage>   resume after a failure; earlier stages are assumed done
#                    (all link calls up to that point still run)
#   --force <stage>  rerun one stage even if it looks complete
#   --dry-run        print the plan and what would be submitted; submit nothing
#   --status         print per-stage completion and exit
#   --no-cancel      do not scancel sibling jobs when a stage fails
#   --help

set -uo pipefail
shopt -s nullglob globstar

BASE_DIR="/sibcb1/jizhelab1/hukaining/202608/OE_WGS"
LINKER="$BASE_DIR/link_stage_outputs.sh"
TUMOR_RAW="$BASE_DIR/WGS"
BLOOD_RAW="$BASE_DIR/WGS_blood/DZOE2026070157-b1/Raw_data"

# ---------------------------------------------------------------------------
# Stage table: id | rel-dir | batcher | done-glob | group | description
# ---------------------------------------------------------------------------
declare -A STAGE_DIR STAGE_BATCHER STAGE_GLOB STAGE_GROUP STAGE_DESC STAGE_ADVISORY

add_stage() {
    STAGE_DIR[$1]="$2"; STAGE_BATCHER[$1]="$3"; STAGE_GLOB[$1]="$4"
    STAGE_GROUP[$1]="$5"; STAGE_DESC[$1]="$6"
}

add_stage 1a    "WGS/1.QC"              batfastp.sh                     "*.R1.clean.fq.gz"                                          G1 "QC tumour (fastp)"
add_stage 1b    "WGS_blood/1.QC"        batfastp.sh                     "*.R1.clean.fq.gz"                                          G1 "QC blood (fastp)"
add_stage 2a    "WGS/2.BWA2_res"        batBWA_MEM2_Picard_markdup.sh   "*.markdup.bam"                                             G2 "BWA-MEM2 + markdup tumour"
add_stage 2b    "WGS_blood/2.BWA2_res"  batBWA_MEM2_Picard_markdup.sh   "*.markdup.bam"                                             G2 "BWA-MEM2 + markdup blood"
add_stage 3a    "WGS/3.BQSR"            batBQSR.sh                      "*.markdup.recal.bam"                                       G3 "BQSR tumour"
add_stage 3b    "WGS_blood/3.BQSR"      batBQSR.sh                      "*.markdup.recal.bam"                                       G3 "BQSR blood"
add_stage 3.1a  "WGS/3.1Pileup_Conta"   batPileup_Conta.sh              "*.pileups.table"                                           G4 "GetPileupSummaries tumour"
add_stage 3.1b  "WGS_blood/3.1Pileup"   batPileup.sh                    "*.pileups.table"                                           G4 "GetPileupSummaries blood"
add_stage 3.2   "3.2CalCont"            batCalCont.sh                   "*.contamination.table"                                     G5 "CalculateContamination"
add_stage 4     "4.SNV"                 batMutect2_BQSR.sh              "*.recal.somatic.SNV.unfiltered.vcf.gz"                     G6 "Mutect2"
add_stage 4_1   "4_1.St_SNV"            batStelka_BQSR.sh               "*st_somatic_BQSR/results/variants/somatic.snvs.vcf.gz"     G6 "Strelka2"
add_stage 6     "6.CNV"                 batCNV_FACETS.sh                "*.facet_out.vcf.gz"                                        G6 "FACETS CNV"
add_stage 7     "7.HLA_type"            batHLAtp.sh                     "**/hla/R1_bestguess_G.txt"                                 G6 "HLA-LA typing"
add_stage 4.f1r2 "4.SNV"                batF1R2Model.sh                 "*.read-orientation-model.tar.gz"                           G7 "F1R2 read-orientation model"
add_stage 5     "5.FilterSNV"           batFilter.sh                    "*.mutect2_strelka.both.PASS.snpeff.missense_frameshift.vcf" G8 "FilterMutectCalls + snpEff"
add_stage 6.1   "6.1phyloWGS"           batPhyloWGS.sh                  "*.top5_tree_ids.tsv"                                       G9 "phyloWGS clone evolution"
add_stage 8     "8.VacQ"                batVacQ.sh                      "quality_output/run_manifest.json"                          G10 "VEP + pVACseq + quality"

# Stage 8 legitimately exits 0 with no output when the both-PASS VCF is empty,
# so its done-glob is advisory rather than required.
STAGE_ADVISORY[8]=1

ORDERED_IDS=(1a 1b 2a 2b 3a 3b 3.1a 3.1b 3.2 4 4_1 6 7 4.f1r2 5 6.1 8)

declare -A GROUP_MEMBERS
for _id in "${ORDERED_IDS[@]}"; do
    GROUP_MEMBERS[${STAGE_GROUP[$_id]}]+="$_id "
done

# Ordered execution plan: groups run in parallel, links wire the next stage.
EXEC_PLAN=(
    "group G1"  "link 1"
    "group G2"  "link 2"
    "group G3"  "link 3"
    "group G4"  "link 3.1"
    "group G5"  "link 3.2"
    "group G6"
    "group G7"  "link 4"  "link 4.f1r2"
    "group G8"  "link 6.1"
    "group G9"  "link 8"
    "group G10"
)

# ---------------------------------------------------------------------------
# Usage / argument parsing
# ---------------------------------------------------------------------------
usage() { awk 'NR>1 && /^#/ {sub(/^# ?/,""); print; next} NR>1 {exit}' "$0"; exit "${1:-0}"; }

SAMPLE=""; FROM=""; FORCE=""; DRY_RUN=false; STATUS=false; DO_CANCEL=true

while (( $# )); do
    case "$1" in
        --from)     FROM="${2:-}"; shift 2 ;;
        --force)    FORCE="${2:-}"; shift 2 ;;
        --dry-run)  DRY_RUN=true; shift ;;
        --status)   STATUS=true; shift ;;
        --no-cancel) DO_CANCEL=false; shift ;;
        --help|-h)  usage 0 ;;
        -*)         echo "Unknown option: $1" >&2; usage 1 ;;
        *)          [[ -z "$SAMPLE" ]] || { echo "Only one SAMPLE_ID allowed" >&2; usage 1; }
                    SAMPLE="$1"; shift ;;
    esac
done

[[ -n "$SAMPLE" ]] || { echo "ERROR: SAMPLE_ID is required" >&2; usage 1; }

P="$BASE_DIR/$SAMPLE"
LOGDIR="$BASE_DIR/logs/$SAMPLE"
JOBDIR="$LOGDIR/jobs"
STATE="$LOGDIR/state.tsv"
LOGFILE="$LOGDIR/run.log"

log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
warn() { log "WARN: $*"; }
die()  { log "ERROR: $*"; log "Resume with: bash $0 $SAMPLE --from <stage>"; exit 1; }

# ---------------------------------------------------------------------------
# Resumability helpers
# ---------------------------------------------------------------------------
record_state() { printf '%s\t%s\t%s\n' "$1" "$2" "$(date '+%F %T')" >> "$STATE"; }

state_rc() {  # <id> -> prints last recorded rc, or nothing
    [[ -f "$STATE" ]] || return 0
    awk -F'\t' -v s="$1" '$1==s{rc=$2} END{if(rc!="")print rc}' "$STATE"
}

# Stage ids strictly before $1 in ORDERED_IDS are skipped when --from is used.
skip_before() {
    [[ -n "$FROM" ]] || return 1
    local id
    for id in "${ORDERED_IDS[@]}"; do
        [[ "$id" == "$FROM" ]] && return 1
        [[ "$id" == "$1" ]] && return 0
    done
    return 1
}

glob_exists() {  # <dir> <pattern>
    local -a m=( "$1"/$2 )
    (( ${#m[@]} > 0 ))
}

stage_done() {  # <id>
    local id="$1" rc
    [[ -n "$FORCE" && "$FORCE" == "$id" ]] && return 1
    skip_before "$id" && return 0
    rc="$(state_rc "$id")"
    [[ "$rc" == "0" ]] && return 0
    [[ -n "$rc" ]] && return 1          # recorded failure -> not done
    glob_exists "$P/${STAGE_DIR[$id]}" "${STAGE_GLOB[$id]}"
}

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
preflight() {
    local missing=() f
    for f in "$TUMOR_RAW/PC_$SAMPLE/PC_$SAMPLE.R1.fastq.gz" \
             "$TUMOR_RAW/PC_$SAMPLE/PC_$SAMPLE.R2.fastq.gz" \
             "$BLOOD_RAW/$SAMPLE/$SAMPLE.R1.fastq.gz" \
             "$BLOOD_RAW/$SAMPLE/$SAMPLE.R2.fastq.gz"; do
        [[ -s "$f" ]] || missing+=("$f")
    done
    if (( ${#missing[@]} )); then
        log "Missing or empty raw FASTQ:"
        printf '  %s\n' "${missing[@]}" | tee -a "$LOGFILE"
        return 1
    fi
    log "preflight: all 4 raw FASTQ present"
    if command -v sinfo >/dev/null 2>&1; then
        sinfo -h -o "%P" >/dev/null 2>&1 || warn "sinfo did not respond"
    fi
    return 0
}

# ---------------------------------------------------------------------------
# Batcher generation + submission
# ---------------------------------------------------------------------------
# Run the batcher WITHOUT `| bash` to only generate .sbatch files, then parse
# the `sbatch <file>` lines it printed. Capturing this run's stdout (rather than
# globbing *.sbatch) is what stops a stale file from an earlier run being
# resubmitted. Never derive sbatch names from a convention: 5.FilterSNV writes
# <ID>.Flter.sbatch (typo) while its job name is <ID>.Filter.
submit_batcher() {  # <id> -> prints generated sbatch basenames; returns 1 on failure
    local id="$1" dir="${STAGE_DIR[$1]}" batcher="${STAGE_BATCHER[$1]}" out line
    [[ -f "$P/$dir/$batcher" ]] || { warn "batcher not found: $P/$dir/$batcher"; return 1; }
    out=$( cd "$P/$dir" && bash "$batcher" ) || { warn "batcher failed: $dir/$batcher"; return 1; }
    while IFS= read -r line; do
        [[ "$line" == "sbatch "* ]] && printf '%s\n' "${line#sbatch }"
    done <<< "$out"
}

jobname_of() {  # <id> <sbatch-file> -> the SLURM job name (and thus .err basename)
    sed -n 's/^#SBATCH --job-name=//p' "$P/${STAGE_DIR[$1]}/$2" | head -1
}

report_failure() {  # <id> <sbatch-basename> <rc>
    local id="$1" f="$2" rc="$3" dir="${STAGE_DIR[$1]}"
    local jn err
    jn="$(jobname_of "$id" "$f")"
    err="$P/$dir/${jn:-unknown}.err"
    echo
    log "==================== STAGE FAILED ===================="
    log "stage     : $id (${STAGE_DESC[$id]})"
    log "directory : $P/$dir"
    log "sbatch    : $dir/$f"
    log "job name  : ${jn:-<unreadable>}"
    log "job id    : $(cat "$JOBDIR/${id}__$(basename "$f").id" 2>/dev/null || echo unknown)"
    log "exit code : $rc"
    log "error log : $err"
    if [[ -f "$err" ]]; then
        if grep -qiE "DUE TO TIME LIMIT" "$err"; then
            log "DIAGNOSIS : WALL-TIME limit exceeded (partition max 72h)"
        elif grep -qiE "oom-kill|Out Of Memory|oom_kill" "$err"; then
            log "DIAGNOSIS : OUT OF MEMORY kill"
        fi
        echo "----- tail -20 $err -----"
        tail -20 "$err"
        echo "-------------------------"
    else
        warn "no .err file at that path; check ${jn:-job}.out in the same directory"
    fi
    log "======================================================"
}

# Run one group: submit every member's jobs concurrently, block until all
# finish, then surface any non-zero exit code.
run_group() {  # <group>
    local g="$1" id f dir rc
    local -a pids=()
    local manifest="$JOBDIR/$g.files"
    : > "$manifest"

    for id in ${GROUP_MEMBERS[$g]}; do
        if stage_done "$id"; then
            log "  [$g] skip $id (${STAGE_DESC[$id]}) -- already complete"
            continue
        fi
        dir="${STAGE_DIR[$id]}"
        # A dry run for a brand-new sample has no batchers yet (scaffolding was
        # intentionally skipped), so report the plan instead of failing.
        if $DRY_RUN && [[ ! -f "$P/$dir/${STAGE_BATCHER[$id]}" ]]; then
            log "  [$g] WOULD SCAFFOLD then SUBMIT $dir (${STAGE_BATCHER[$id]})   (${STAGE_DESC[$id]})"
            continue
        fi
        local -a files=()
        local sbout
        sbout="$(submit_batcher "$id")" \
            || die "stage $id (${STAGE_DESC[$id]}): could not generate sbatch files in $P/$dir"
        [[ -n "$sbout" ]] && mapfile -t files <<< "$sbout"
        if (( ${#files[@]} == 0 )); then
            die "stage $id (${STAGE_DESC[$id]}): batcher generated no jobs. Its input glob matched nothing in $P/$dir -- are the upstream symlinks wired?"
        fi
        for f in "${files[@]}"; do
            [[ -f "$P/$dir/$f" ]] || die "stage $id: batcher announced a missing file: $dir/$f"
            printf '%s\t%s\t%s\n' "$id" "$dir" "$f" >> "$manifest"
            if $DRY_RUN; then
                log "  [$g] WOULD SUBMIT $dir/$f   (${STAGE_DESC[$id]})"
                continue
            fi
            log "  [$g] submit $dir/$f   (${STAGE_DESC[$id]})"
            # The generated sbatch uses RELATIVE paths (input files, --output,
            # --error), so SLURM must start the job in the stage directory --
            # that is what `cd <stage> && bash batcher | bash` gives for free.
            ( cd "$P/$dir" || exit 1
              out=$(sbatch --wait --parsable "$f" 2>&1 ); r=$?
              printf '%s\n' "${out%%;*}" > "$JOBDIR/${id}__$(basename "$f").id"
              printf '%s\n' "$r"        > "$JOBDIR/${id}__$(basename "$f").rc"
              printf '%s\n' "$out"       > "$JOBDIR/${id}__$(basename "$f").sub" ) &
            pids+=("$!")
        done
    done

    $DRY_RUN && return 0
    (( ${#pids[@]} )) || { log "  [$g] nothing to submit"; return 0; }

    local pid
    for pid in "${pids[@]}"; do wait "$pid" || true; done

    # Collect outcomes; report the first failure with full diagnostics.
    local failed_id="" failed_file="" failed_rc=0
    while IFS=$'\t' read -r id dir f; do
        rc="$(cat "$JOBDIR/${id}__$(basename "$f").rc" 2>/dev/null || echo 1)"
        if (( rc != 0 )); then
            log "  [$g] FAILED $dir/$f (exit $rc)"
            [[ -z "$failed_id" ]] && { failed_id="$id"; failed_file="$f"; failed_rc="$rc"; }
        else
            log "  [$g] ok     $dir/$f"
        fi
    done < "$manifest"

    if [[ -n "$failed_id" ]]; then
        if $DO_CANCEL; then
            while IFS=$'\t' read -r id dir f; do
                [[ "$id" == "$failed_id" && "$f" == "$failed_file" ]] && continue
                local jid; jid="$(cat "$JOBDIR/${id}__$(basename "$f").id" 2>/dev/null)"
                [[ -n "$jid" ]] && scancel "$jid" 2>/dev/null && log "  cancelled sibling job $jid"
            done < "$manifest"
        fi
        # Record the failure BEFORE dying: a failed stage may have written a
        # partial output, and on resume the done-glob alone would call it done.
        # stage_done() treats any recorded non-zero rc as "not done".
        record_state "$failed_id" "$failed_rc"
        report_failure "$failed_id" "$failed_file" "$failed_rc"
        die "stage $failed_id failed; chain stopped"
    fi

    # Exit codes are not sufficient on their own: no generated sbatch sets
    # `set -e`, so rc reflects only the last command. Require the output too.
    while IFS=$'\t' read -r id dir f; do
        if ! glob_exists "$P/$dir" "${STAGE_GLOB[$id]}"; then
            if [[ -n "${STAGE_ADVISORY[$id]:-}" ]]; then
                warn "stage $id produced no ${STAGE_GLOB[$id]} -- ok if the both-PASS VCF was empty"
            else
                # Sentinel rc 99 so a resume retries this rather than calling
                # it done off a stale or partial file.
                record_state "$id" 99
                report_failure "$id" "$f" 0
                die "stage $id exited 0 but produced no ${STAGE_GLOB[$id]} in $P/$dir"
            fi
        fi
        record_state "$id" 0
    done < "$manifest"
}

run_link() {  # <stage>
    local stage="$1"
    if $DRY_RUN; then log "  link $stage (would run link_stage_outputs.sh $SAMPLE $stage)"; return 0; fi
    log "  link $stage"
    bash "$LINKER" "$SAMPLE" "$stage" >>"$LOGFILE" 2>&1 \
        || die "link_stage_outputs.sh $SAMPLE $stage failed"
}

# ---------------------------------------------------------------------------
# Status / plan
# ---------------------------------------------------------------------------
print_status() {
    printf '%-8s %-34s %-10s %s\n' STAGE DESCRIPTION STATUS GLOB
    local id rc mark
    for id in "${ORDERED_IDS[@]}"; do
        rc="$(state_rc "$id")"
        if [[ "$rc" == "0" ]]; then mark="ok"
        elif [[ -n "$rc" ]]; then mark="FAILED($rc)"
        elif glob_exists "$P/${STAGE_DIR[$id]}" "${STAGE_GLOB[$id]}"; then mark="output-present"
        else mark="-"
        fi
        printf '%-8s %-34s %-10s %s\n' "$id" "${STAGE_DESC[$id]}" "$mark" "${STAGE_GLOB[$id]}"
    done
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
# --status is a read-only query: answer it without creating log files.
if $STATUS; then
    [[ -d "$P" ]] || { echo "ERROR: no such patient directory: $P" >&2; exit 1; }
    print_status
    exit 0
fi

# Validate arguments before creating anything.
if [[ -n "$FROM" && -z "${STAGE_DIR[$FROM]:-}" ]]; then
    echo "ERROR: --from '$FROM' is not a known stage. Valid: ${ORDERED_IDS[*]}" >&2
    exit 1
fi
if [[ -n "$FORCE" && -z "${STAGE_DIR[$FORCE]:-}" ]]; then
    echo "ERROR: --force '$FORCE' is not a known stage. Valid: ${ORDERED_IDS[*]}" >&2
    exit 1
fi

mkdir -p "$JOBDIR"
# Mirror everything to the log from here on.
exec > >(tee -a "$LOGFILE") 2>&1

log "=== run_pipeline.sh $SAMPLE ==="
log "host=$(hostname)  pid=$$  tmux=${TMUX:+yes}${TMUX:-no}"
if [[ -z "${TMUX:-}" ]]; then
    warn "not running inside tmux -- if this session dies the driver stops."
    warn "strongly recommended: tmux new -s $SAMPLE"
fi
$DO_CANCEL || log "note: --no-cancel set; sibling jobs will be left running on failure"

preflight || die "raw FASTQ preflight failed"

if [[ -d "$P" && -z "$FROM" && -f "$STATE" ]]; then
    warn "$P already has run state. Resuming; stages recorded ok will be skipped."
fi

if ! $DRY_RUN; then
    running=$(squeue -h -u "$USER" -o "%j" 2>/dev/null | grep -c "$SAMPLE" || true)
    if (( running > 0 )); then
        warn "$running queued/running SLURM job(s) mention '$SAMPLE'."
        warn "If a previous driver is still alive, do not start a second one."
    fi

    log "--- scaffolding ---"
    bash "$BASE_DIR/setup_pipeline.sh"    "$SAMPLE" >>"$LOGFILE" 2>&1 \
        || die "setup_pipeline.sh $SAMPLE failed"
    bash "$BASE_DIR/setup_rest_steps.sh"  "$SAMPLE" >>"$LOGFILE" 2>&1 \
        || die "setup_rest_steps.sh $SAMPLE failed"
    log "scaffold complete"
else
    log "--- dry run: skipping scaffolding ---"
fi

for step in "${EXEC_PLAN[@]}"; do
    read -r kind arg <<< "$step"
    case "$kind" in
        group) run_group "$arg" ;;
        link)  run_link  "$arg" ;;
    esac
done

log "=== DONE: $SAMPLE completed all stages ==="
