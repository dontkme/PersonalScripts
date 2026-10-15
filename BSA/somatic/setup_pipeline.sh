#!/bin/bash
# Setup pipeline directories, FASTQ symlinks, and batcher scripts for all
# remaining OE_WGS samples (skip BE0728, E88052, E88354 which already exist).
#
# Usage:
#   bash setup_pipeline.sh           # Create dirs, symlinks, and batcher scripts
#   bash setup_pipeline.sh --submit  # Also auto-submit all sbatch jobs
#
# The script creates per-sample pipeline dirs with WGS/ and WGS_blood/
# subdirectories, each containing stage dirs with batcher scripts.
# Run each batcher from its stage directory to generate sbatch files:
#   cd B1C391/WGS/1.QC && bash batfastp.sh

set -euo pipefail

BASE_DIR="/sibcb1/jizhelab1/hukaining/202608/OE_WGS"
TUMOR_RAW="${BASE_DIR}/WGS"
BLOOD_RAW="${BASE_DIR}/WGS_blood/DZOE2026070157-b1/Raw_data"

# Tool paths
PICARD="/sibcb1/jizhelab1/hukaining/Biosoft/picard.jar"
GATK="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
REF_FA="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
REF_BWA2="/sibcb1/jizhelab1/hukaining/genome/hg38.BWA_MEM2_index/GRCh38.primary_assembly.genome.fa"
DBSNP="/sibcb1/jizhelab1/hukaining/genome/dbsnp/b151/hg38p7/00-All.vcf.gz"
MILLS_INDELS="/sibcb1/jizhelab1/hukaining/genome/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# Samples that don't have pipeline dirs yet. Override by passing sample IDs
# as positional args: `bash setup_pipeline.sh NEWSAMPLE`
SAMPLES=(B1C391 BDZ717 BEC478 BEF141 BEG091 BEG323 BEH047 E89007)

AUTO_SUBMIT=false
ARGS=()
for a in "$@"; do
    if [[ "$a" == "--submit" ]]; then
        AUTO_SUBMIT=true
    else
        ARGS+=("$a")
    fi
done
(( ${#ARGS[@]} )) && SAMPLES=("${ARGS[@]}")

###############################################################################
# Helper: write a batcher script
###############################################################################
write_batcher() {
    local dest="$1"
    local content="$2"
    cat > "$dest" <<'BATCHER_EOF'
BATCHER_EOF
    # Use printf to inject the content safely
    printf '%s\n' "$content" >> "$dest"
    chmod +x "$dest"
}

###############################################################################
# 1.QC: fastp batcher
###############################################################################
gen_fastp_batcher() {
    local sample="$1"          # e.g., B1C391
    local sample_type="$2"     # Tumor or Blood
    local prefix
    if [[ "$sample_type" == "Tumor" ]]; then
        prefix="PC_${sample}"
    else
        prefix="${sample}"
    fi

    cat <<EOF
#!/bin/bash
for file in \$(ls *R1.fastq.gz); do
    prefix=\$(echo "\$file" | sed 's/\(.*\)\\.R1\\.fastq\\.gz/\\1/g')
    file1="\${prefix}.R1.fastq.gz"
    file2="\${prefix}.R2.fastq.gz"
    short="\${prefix}.${sample_type}"
    echo "#!/bin/bash" >"\$short.sbatch"
    echo "#SBATCH --job-name=\$short.fastp" >>"\$short.sbatch"
    echo "#SBATCH --output=\$short.fastp.out" >>"\$short.sbatch"
    echo "#SBATCH --error=\$short.fastp.err" >>"\$short.sbatch"
    echo "#SBATCH --partition=default1" >>"\$short.sbatch"
    echo "#SBATCH --ntasks=1" >>"\$short.sbatch"
    echo "#SBATCH --cpus-per-task=4" >>"\$short.sbatch"
    echo "#SBATCH --mem=64G" >>"\$short.sbatch"
    echo "fastp -w 4 -i \$file1 -I \$file2 -o \$short.R1.clean.fq.gz -O \$short.R2.clean.fq.gz --detect_adapter_for_pe --n_base_limit 4 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --length_required 75 --average_qual 15 -h \$short.fastp_report.html -j \$short.fastp_report.json" >>"\$short.sbatch"
    echo "sbatch \$short.sbatch"
done
EOF
}

###############################################################################
# 2.BWA2_res: BWA-MEM2 + Picard MarkDuplicates batcher
###############################################################################
gen_bwa2_batcher() {
    cat <<'EOF'
#!/bin/bash
picardpath="/sibcb1/jizhelab1/hukaining/Biosoft/picard.jar"
refbwa="/sibcb1/jizhelab1/hukaining/genome/hg38.BWA_MEM2_index/GRCh38.primary_assembly.genome.fa"
bwancpu=12
samtoolsncpu=4
for i in $(ls *R1.clean.fq.gz); do
    short=$(echo $i | sed 's/\(.*\)\.R1\.clean\.fq\.gz/\1/')
    file1="$short.R1.clean.fq.gz"
    file2="$short.R2.clean.fq.gz"
    jobname="$short.BWA.mem"
    sbatchname="$short.BWA.mem.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=default1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=16" >>$sbatchname
    echo "#SBATCH --mem=120G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    printf '%s\n' "bwa-mem2 mem -R '@RG\\tID:${short}\\tSM:${short}\\tLB:${short}\\tPL:ILLUMINA' -t ${bwancpu} -K 100000000 \"${refbwa}\" \"${file1}\" \"${file2}\" | samtools sort -@ ${samtoolsncpu} -m 4G -o \"${short}.sorted.bam\" -" >> "$sbatchname"
    echo "samtools index $short.sorted.bam" >>$sbatchname
    echo "samtools flagstat $short.sorted.bam >$short.flagstat" >>$sbatchname
    echo "java -Xmx110G -jar $picardpath MarkDuplicates I=$short.sorted.bam O=$short.markdup.bam M=$short.rmdup.metrics.txt REMOVE_DUPLICATES=false CREATE_INDEX=true" >>$sbatchname
    echo "samtools flagstat $short.markdup.bam >$short.markdup.flagstat" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 3.BQSR: GATK BaseRecalibrator + ApplyBQSR batcher
###############################################################################
gen_bqsr_batcher() {
    cat <<'EOF'
#!/bin/bash
gatkpath="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
dbsnppath="/sibcb1/jizhelab1/hukaining/genome/dbsnp/b151/hg38p7/00-All.vcf.gz"
milles_indelspath="/sibcb1/jizhelab1/hukaining/genome/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
reffa="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
for i in $(ls *.markdup.bam); do
    short=$(echo $i | sed 's/\(.*\)\.markdup\.bam/\1/')
    file=$i
    jobname="$short.BQSR"
    sbatchname="$short.BQSR.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=default1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=12" >>$sbatchname
    echo "#SBATCH --mem=64G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    outf="$short.markdup.recal"
    echo "$gatkpath --java-options \"-Xmx64g -XX:+UseG1GC -XX:ParallelGCThreads=4 -XX:ConcGCThreads=2\" BaseRecalibrator -R $reffa --known-sites $dbsnppath --known-sites $milles_indelspath -I $file -O $outf.table" >>$sbatchname
    echo "$gatkpath --java-options \"-Xmx64g -XX:+UseG1GC -XX:ParallelGCThreads=4 -XX:ConcGCThreads=2\" ApplyBQSR -R $reffa -I $file -O $outf.bam --bqsr-recal-file $outf.table" >>$sbatchname
    echo "samtools flagstat $outf.bam >$outf.bam.flagstat" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# Main
###############################################################################
echo "=== OE_WGS Pipeline Setup ==="
echo "Samples: ${SAMPLES[*]}"
echo "Auto-submit: $AUTO_SUBMIT"
echo ""

for sample in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${BASE_DIR}/${sample}"

    echo "--- Setting up ${sample} ---"

    # ---- Create directory structure ----
    mkdir -p "${SAMPLE_DIR}/WGS/1.QC"
    mkdir -p "${SAMPLE_DIR}/WGS/2.BWA2_res"
    mkdir -p "${SAMPLE_DIR}/WGS/3.BQSR"
    mkdir -p "${SAMPLE_DIR}/WGS_blood/1.QC"
    mkdir -p "${SAMPLE_DIR}/WGS_blood/2.BWA2_res"
    mkdir -p "${SAMPLE_DIR}/WGS_blood/3.BQSR"

    # ---- Symlink tumor FASTQ into WGS/ parent (stage dirs point here via ../) ----
    ln -sf "../../WGS/PC_${sample}/PC_${sample}.R1.fastq.gz" "${SAMPLE_DIR}/WGS/PC_${sample}.R1.fastq.gz"
    ln -sf "../../WGS/PC_${sample}/PC_${sample}.R2.fastq.gz" "${SAMPLE_DIR}/WGS/PC_${sample}.R2.fastq.gz"

    # ---- Symlink blood FASTQ into WGS_blood/ parent ----
    ln -sf "../../WGS_blood/DZOE2026070157-b1/Raw_data/${sample}/${sample}.R1.fastq.gz" "${SAMPLE_DIR}/WGS_blood/${sample}.R1.fastq.gz"
    ln -sf "../../WGS_blood/DZOE2026070157-b1/Raw_data/${sample}/${sample}.R2.fastq.gz" "${SAMPLE_DIR}/WGS_blood/${sample}.R2.fastq.gz"

    # ---- Symlink raw FASTQ into step 1 stage dirs only ----
    # Step 2 and 3 don't need FASTQ — they consume the previous step's output
    # (use link_stage_outputs.sh after each step completes to link outputs forward)
    ln -sf "../PC_${sample}.R1.fastq.gz" "${SAMPLE_DIR}/WGS/1.QC/PC_${sample}.R1.fastq.gz"
    ln -sf "../PC_${sample}.R2.fastq.gz" "${SAMPLE_DIR}/WGS/1.QC/PC_${sample}.R2.fastq.gz"
    ln -sf "../${sample}.R1.fastq.gz" "${SAMPLE_DIR}/WGS_blood/1.QC/${sample}.R1.fastq.gz"
    ln -sf "../${sample}.R2.fastq.gz" "${SAMPLE_DIR}/WGS_blood/1.QC/${sample}.R2.fastq.gz"

    # ---- Write batcher scripts ----
    gen_fastp_batcher "$sample" "Tumor"  > "${SAMPLE_DIR}/WGS/1.QC/batfastp.sh"
    gen_fastp_batcher "$sample" "Blood"  > "${SAMPLE_DIR}/WGS_blood/1.QC/batfastp.sh"
    gen_bwa2_batcher                     > "${SAMPLE_DIR}/WGS/2.BWA2_res/batBWA_MEM2_Picard_markdup.sh"
    gen_bwa2_batcher                     > "${SAMPLE_DIR}/WGS_blood/2.BWA2_res/batBWA_MEM2_Picard_markdup.sh"
    gen_bqsr_batcher                     > "${SAMPLE_DIR}/WGS/3.BQSR/batBQSR.sh"
    gen_bqsr_batcher                     > "${SAMPLE_DIR}/WGS_blood/3.BQSR/batBQSR.sh"

    chmod +x "${SAMPLE_DIR}/WGS/1.QC/batfastp.sh"
    chmod +x "${SAMPLE_DIR}/WGS_blood/1.QC/batfastp.sh"
    chmod +x "${SAMPLE_DIR}/WGS/2.BWA2_res/batBWA_MEM2_Picard_markdup.sh"
    chmod +x "${SAMPLE_DIR}/WGS_blood/2.BWA2_res/batBWA_MEM2_Picard_markdup.sh"
    chmod +x "${SAMPLE_DIR}/WGS/3.BQSR/batBQSR.sh"
    chmod +x "${SAMPLE_DIR}/WGS_blood/3.BQSR/batBQSR.sh"

    echo "  Created: ${SAMPLE_DIR}/WGS/{1.QC,2.BWA2_res,3.BQSR}"
    echo "  Created: ${SAMPLE_DIR}/WGS_blood/{1.QC,2.BWA2_res,3.BQSR}"
    echo "  Symlinks: tumor FASTQ -> WGS/PC_${sample}.R{1,2}.fastq.gz (and into 1.QC/)"
    echo "  Symlinks: blood FASTQ -> WGS_blood/${sample}.R{1,2}.fastq.gz (and into 1.QC/)"
    echo "  NOTE: After step 1 completes, run link_stage_outputs.sh to symlink clean.fq.gz into step 2."
    echo "  NOTE: After step 2 completes, run link_stage_outputs.sh to symlink markdup.bam into step 3."

    # ---- Optionally submit jobs ----
    if $AUTO_SUBMIT; then
        echo "  Submitting jobs..."
        (cd "${SAMPLE_DIR}/WGS/1.QC" && bash batfastp.sh | bash)
        (cd "${SAMPLE_DIR}/WGS_blood/1.QC" && bash batfastp.sh | bash)
        echo "  Submitted fastp jobs for ${sample} (tumor + blood)"
        echo "  NOTE: Run BWA2 after fastp completes, then BQSR after BWA2 completes."
    fi

    echo ""
done

echo "=== Done ==="
echo "Pipeline directories created for ${#SAMPLES[@]} samples."
echo ""
echo "Next steps (per sample, in order):"
echo "  # Step 1: QC"
echo "  cd ${BASE_DIR}/B1C391/WGS/1.QC && bash batfastp.sh"
echo "  cd ${BASE_DIR}/B1C391/WGS_blood/1.QC && bash batfastp.sh"
echo ""
echo "  # After step 1 finishes, link outputs forward:"
echo "  bash link_stage_outputs.sh B1C391 1   # symlinks clean.fq.gz into 2.BWA2_res/"
echo ""
echo "  # Step 2: Alignment"
echo "  cd ${BASE_DIR}/B1C391/WGS/2.BWA2_res && bash batBWA_MEM2_Picard_markdup.sh"
echo "  cd ${BASE_DIR}/B1C391/WGS_blood/2.BWA2_res && bash batBWA_MEM2_Picard_markdup.sh"
echo ""
echo "  # After step 2 finishes, link outputs forward:"
echo "  bash link_stage_outputs.sh B1C391 2   # symlinks markdup.bam into 3.BQSR/"
echo ""
echo "  # Step 3: BQSR"
echo "  cd ${BASE_DIR}/B1C391/WGS/3.BQSR && bash batBQSR.sh"
echo "  cd ${BASE_DIR}/B1C391/WGS_blood/3.BQSR && bash batBQSR.sh"
