#!/bin/bash
# Setup downstream pipeline stages (3.1-8) for all OE_WGS samples.
# Depends on stages 1-3 already being set up by setup_pipeline.sh.
#
# Usage:
#   bash setup_rest_steps.sh           # Create dirs and batcher scripts
#
# After each stage completes, use link_stage_outputs.sh to wire inputs forward.

set -euo pipefail

BASE_DIR="/sibcb1/jizhelab1/hukaining/202608/OE_WGS"
#ALL_SAMPLES=(B1C391 BDZ717 BEC478 BEF141 BEG091 BEG323 BEH047 E89007)
ALL_SAMPLES=(BE0728 E88354 E88052 B1C391 BDZ717 BEC478 BEF141 BEG091 BEG323 BEH047 E89007)

# Override by passing sample IDs as positional args:
# `bash setup_rest_steps.sh NEWSAMPLE`
(( $# )) && ALL_SAMPLES=("$@")

###############################################################################
# 3.1Pileup_Conta (WGS/tumor side): GATK GetPileupSummaries on recal.bam
###############################################################################
gen_pileup_tumor_batcher() {
    cat <<'EOF'
#!/bin/bash
gatkpath="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
AF_gnomadbipath="/sibcb1/jizhelab1/hukaining/genome/gnomeAD/af-only-gnomad.hg38.main.common-biallelic.vcf.gz"
reffa="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
ncpu=12
for i in $(ls *.recal.bam); do
    short=$(echo $i | sed 's/\(.*\)\.markdup\.recal\.bam/\1/')
    file=$i
    jobname="$short.Pileup"
    sbatchname="$short.Pileup.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=default1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=120G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    outfp=$short.pileups
    echo "$gatkpath --java-options \"-Xmx110g -XX:+UseG1GC -XX:ParallelGCThreads=12 -XX:ConcGCThreads=12\" GetPileupSummaries -R $reffa -I $file -V $AF_gnomadbipath -L $AF_gnomadbipath -O $outfp.table" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 3.1Pileup (WGS_blood/blood side): GATK GetPileupSummaries on recal.bam
###############################################################################
gen_pileup_blood_batcher() {
    cat <<'EOF'
#!/bin/bash
gatkpath="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
AF_gnomadbipath="/sibcb1/jizhelab1/hukaining/genome/gnomeAD/af-only-gnomad.hg38.main.common-biallelic.vcf.gz"
reffa="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
ncpu=12
for i in $(ls *.recal.bam); do
    short=$(echo $i | sed 's/\(.*\)\.markdup\.recal\.bam/\1/')
    file=$i
    jobname="$short.Pileup"
    sbatchname="$short.Pileup.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=default1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=120G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    outfp=$short.pileups
    echo "$gatkpath --java-options \"-Xmx110g -XX:+UseG1GC -XX:ParallelGCThreads=12 -XX:ConcGCThreads=12\" GetPileupSummaries -R $reffa -I $file -V $AF_gnomadbipath -L $AF_gnomadbipath -O $outfp.table" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 3.2CalCont (top-level): GATK CalculateContamination
###############################################################################
gen_calcont_batcher() {
    cat <<'EOF'
#!/bin/bash
gatkpath="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
ncpu=16
for i in $(ls *.Blood.pileups.table); do
    short=$(echo $i | sed 's/\(.*\)\.Blood\.pileups\.table/\1/')
    fileB=$i
    fileT=PC_$short.Tumor.pileups.table
    jobname="$short.CalCont"
    sbatchname="$short.CalCont.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=gpu1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=64G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    outfc=$short.contamination
    echo "$gatkpath --java-options \"-Xmx64g -XX:+UseG1GC -XX:ParallelGCThreads=$ncpu -XX:ConcGCThreads=$ncpu\" CalculateContamination -I $fileT -matched $fileB -O $outfc.table --tumor-segmentation $short.segments.table" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 4.SNV (top-level): GATK Mutect2
###############################################################################
gen_mutect2_batcher() {
    cat <<'EOF'
#!/bin/bash
gatkpath="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
AF_gnomadpath="/sibcb1/jizhelab1/hukaining/genome/gnomeAD/af-only-gnomad.hg38.vcf.gz"
reffa="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
ncpu=12
for i in $(ls *Blood.markdup.recal.bam); do
    short=$(echo $i | sed 's/\(.*\)\.Blood\.markdup\.recal\.bam/\1/')
    fileT=PC_$short.Tumor.markdup.recal.bam
    fileN=$short.Blood.markdup.recal.bam
    SM=$short.Blood
    jobname="$short.SNV.recal"
    sbatchname="$short.SNV.recal.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=gpu1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=120G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    outf=$short.recal.somatic.SNV
    echo "$gatkpath --java-options \"-Xmx110g -XX:+UseG1GC -XX:ParallelGCThreads=$ncpu -XX:ConcGCThreads=$ncpu\" Mutect2 --native-pair-hmm-threads $ncpu -R $reffa -I $fileT -I $fileN -normal $SM --germline-resource $AF_gnomadpath -O $outf.unfiltered.vcf.gz --f1r2-tar-gz $outf.f1r2.tar.gz --bam-output $outf.mutect2.bamout.bam" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 4.SNV F1R2 model (top-level): GATK LearnReadOrientationModel
# Runs AFTER Mutect2 (consumes the f1r2.tar.gz Mutect2 emits).
###############################################################################
gen_f1r2_batcher() {
    cat <<'EOF'
#!/bin/bash
gatkpath="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
ncpu=12
for i in $(ls *Blood.markdup.recal.bam); do
    short=$(echo $i | sed 's/\(.*\)\.Blood\.markdup\.recal\.bam/\1/')
    jobname="$short.SNV.F1R2"
    sbatchname="$short.SNV.F1R2.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=gpu1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=64G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    outf=$short.recal.somatic.SNV
    echo "$gatkpath --java-options \"-Xmx64g -XX:+UseG1GC -XX:ParallelGCThreads=$ncpu -XX:ConcGCThreads=$ncpu\" LearnReadOrientationModel -I $outf.f1r2.tar.gz -O $outf.read-orientation-model.tar.gz" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 4_1.St_SNV (top-level): Strelka2 somatic
###############################################################################
gen_strelka_batcher() {
    cat <<'EOF'
#!/bin/bash
stpath="/sibcb1/jizhelab1/hukaining/Biosoft/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py"
reffa="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
ncpu=12
for i in $(ls *Blood.markdup.recal.bam); do
    short=$(echo $i | sed 's/\(.*\)\.Blood\.markdup\.recal\.bam/\1/')
    fileT=PC_$short.Tumor.markdup.recal.bam
    fileN=$short.Blood.markdup.recal.bam
    jobname="$short.ST_SNV_BQSR"
    sbatchname="$short.ST_SNV_BQSR.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=gpu1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=120G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    echo "mamba activate py27" >>$sbatchname
    outf=$short.st_somatic_BQSR
    echo "$stpath --normalBam $fileN --tumorBam $fileT --referenceFasta $reffa --runDir ${outf}" >>$sbatchname
    echo "${outf}/runWorkflow.py -m local -j $ncpu" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 5.FilterSNV (top-level): FilterMutectCalls + bcftools PASS + Strelka2 merge
#                          + snpEff annotation (Mutect2-only and both-callers)
###############################################################################
gen_filter_batcher() {
    cat <<'EOF'
#!/bin/bash
gatkpath="/sibcb1/jizhelab1/hukaining/Biosoft/gatk-4.7.0.0/gatk"
snpEffpath="/sibcb1/jizhelab1/hukaining/Biosoft/snpEff/snpEff.jar"
SnpSiftpath="/sibcb1/jizhelab1/hukaining/Biosoft/snpEff/SnpSift.jar"
reffa="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
ncpu=16
for i in $(ls *.recal.somatic.SNV.unfiltered.vcf.gz); do
    short=$(echo $i | sed 's/\(.*\)\.recal\.somatic\.SNV\.unfiltered\.vcf\.gz/\1/')
    jobname="$short.Filter"
    sbatchname="$short.Flter.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=gpu1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=64G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    # ---- Mutect2 filtering ----
    echo "$gatkpath --java-options \"-Xmx64g -XX:+UseG1GC -XX:ParallelGCThreads=$ncpu -XX:ConcGCThreads=$ncpu\" FilterMutectCalls -R $reffa -V $i --stats $i.stats --orientation-bias-artifact-priors $short.recal.somatic.SNV.read-orientation-model.tar.gz -contamination-table $short.contamination.table -O $short.mutect2.filtered.vcf.gz --tumor-segmentation $short.segments.table --filtering-stats $short.filtering.stats" >>$sbatchname
    echo "bcftools view -f PASS -Oz -o $short.mutect2.filtered.PASS.vcf.gz $short.mutect2.filtered.vcf.gz" >>$sbatchname
    echo "bcftools norm -f $reffa -m -any -Oz -o $short.mutect2.filtered.PASS.norm.vcf.gz $short.mutect2.filtered.PASS.vcf.gz" >>$sbatchname
    echo "bcftools index -t $short.mutect2.filtered.PASS.norm.vcf.gz" >>$sbatchname
    # ---- Strelka2: merge SNV+indel VCFs, normalize to match Mutect2 ----
    echo "bcftools concat -a -Oz -o $short.strelka.somatic.vcf.gz $short.strelka.somatic.snvs.vcf.gz $short.strelka.somatic.indels.vcf.gz" >>$sbatchname
    echo "bcftools norm -f $reffa -m -any -Oz -o $short.strelka.somatic.norm.vcf.gz $short.strelka.somatic.vcf.gz" >>$sbatchname
    # ---- Tag each Mutect2 record with Strelka2 PASS/NOPASS ----
    echo "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n' $short.strelka.somatic.norm.vcf.gz | awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$4,(\$5==\"PASS\"?\"PASS\":\"NOPASS\")}' | sort -k1,1 -k2,2n | bgzip -c > $short.strelka.status.tab.gz" >>$sbatchname
    echo "tabix -s1 -b2 -e2 $short.strelka.status.tab.gz" >>$sbatchname
    echo "printf '%s\n' '##INFO=<ID=STRELKA,Number=1,Type=String,Description=\"Strelka2 somatic status: PASS or NOPASS; tag absent means not called by Strelka2\">' > $short.strelka.hdr" >>$sbatchname
    echo "bcftools annotate -a $short.strelka.status.tab.gz -h $short.strelka.hdr -c CHROM,POS,REF,ALT,INFO/STRELKA -Oz -o $short.mutect2.filtered.PASS.norm.strelka.vcf.gz $short.mutect2.filtered.PASS.norm.vcf.gz" >>$sbatchname
    echo "bcftools index -t $short.mutect2.filtered.PASS.norm.strelka.vcf.gz" >>$sbatchname
    # ---- VCF requiring both callers to PASS ----
    echo "bcftools view -i 'STRELKA=\"PASS\"' -Oz -o $short.mutect2_strelka.both.PASS.vcf.gz $short.mutect2.filtered.PASS.norm.strelka.vcf.gz" >>$sbatchname
    echo "bcftools index -t $short.mutect2_strelka.both.PASS.vcf.gz" >>$sbatchname
    # ---- snpEff annotation: Mutect2-only ----
    echo "java -Xmx16G -jar $snpEffpath GRCh38.115 -stats $short.mutect2.filtered.PASS.norm.snpeff.html $short.mutect2.filtered.PASS.norm.vcf.gz > $short.mutect2.filtered.PASS.norm.snpeff.vcf" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath filter \"(ANN[*].BIOTYPE = 'protein_coding') & ((ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'frameshift_variant'))\" $short.mutect2.filtered.PASS.norm.snpeff.vcf > $short.mutect2.filtered.PASS.norm.snpeff.missense_frameshift.vcf" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath extractFields -s \",\" -e \".\" $short.mutect2.filtered.PASS.norm.snpeff.missense_frameshift.vcf CHROM POS ID REF ALT QUAL FILTER \"ANN[*].GENE\" \"ANN[*].EFFECT\" \"ANN[*].IMPACT\" \"ANN[*].FEATUREID\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].RANK\" TLOD > $short.mutect2.filtered.PASS.norm.snpeff.missense_frameshift.tsv" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath filter \"(ANN[*].BIOTYPE = 'protein_coding') & ((ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'stop_gained') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'start_lost') | (ANN[*].EFFECT has 'stop_lost'))\" $short.mutect2.filtered.PASS.norm.snpeff.vcf > $short.mutect2.filtered.PASS.norm.snpeff.protein_damaging.vcf" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath extractFields -s \",\" -e \".\" $short.mutect2.filtered.PASS.norm.snpeff.protein_damaging.vcf CHROM POS ID REF ALT QUAL FILTER \"ANN[*].GENE\" \"ANN[*].EFFECT\" \"ANN[*].IMPACT\" \"ANN[*].FEATUREID\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].RANK\" TLOD > $short.mutect2.filtered.PASS.norm.snpeff.protein_damaging.tsv" >>$sbatchname
    # ---- snpEff annotation: Mutect2 + Strelka2 both PASS ----
    echo "java -Xmx16G -jar $snpEffpath GRCh38.115 -stats $short.mutect2_strelka.both.PASS.snpeff.html $short.mutect2_strelka.both.PASS.vcf.gz > $short.mutect2_strelka.both.PASS.snpeff.vcf" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath filter \"(ANN[*].BIOTYPE = 'protein_coding') & ((ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'frameshift_variant'))\" $short.mutect2_strelka.both.PASS.snpeff.vcf > $short.mutect2_strelka.both.PASS.snpeff.missense_frameshift.vcf" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath extractFields -s \",\" -e \".\" $short.mutect2_strelka.both.PASS.snpeff.missense_frameshift.vcf CHROM POS ID REF ALT QUAL FILTER \"ANN[*].GENE\" \"ANN[*].EFFECT\" \"ANN[*].IMPACT\" \"ANN[*].FEATUREID\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].RANK\" TLOD > $short.mutect2_strelka.both.PASS.snpeff.missense_frameshift.tsv" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath filter \"(ANN[*].BIOTYPE = 'protein_coding') & ((ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'stop_gained') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'start_lost') | (ANN[*].EFFECT has 'stop_lost'))\" $short.mutect2_strelka.both.PASS.snpeff.vcf > $short.mutect2_strelka.both.PASS.snpeff.protein_damaging.vcf" >>$sbatchname
    echo "java -Xmx16G -jar $SnpSiftpath extractFields -s \",\" -e \".\" $short.mutect2_strelka.both.PASS.snpeff.protein_damaging.vcf CHROM POS ID REF ALT QUAL FILTER \"ANN[*].GENE\" \"ANN[*].EFFECT\" \"ANN[*].IMPACT\" \"ANN[*].FEATUREID\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].RANK\" TLOD > $short.mutect2_strelka.both.PASS.snpeff.protein_damaging.tsv" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 6.CNV (top-level): FACETS allele-specific copy number + purity/ploidy
###############################################################################
gen_cnv_batcher() {
    cat <<'EOF'
#!/bin/bash
cnvfacets="/sibcb1/jizhelab1/hukaining/miniforge3/envs/py27/bin/cnv_facets.R"
#dbsnppath="/sibcb1/jizhelab1/hukaining/genome/dbsnp/b151/hg38p7/00-All.vcf.gz"
dbsnppath="/sibcb1/jizhelab1/hukaining/genome/dbsnp/00-common_all.chr.vcf.gz"
ncpu=12
for i in $(ls *Blood.markdup.recal.bam); do
    short=$(echo $i | sed 's/\(.*\)\.Blood\.markdup\.recal\.bam/\1/')
    fileT=PC_$short.Tumor.markdup.recal.bam
    fileN=$short.Blood.markdup.recal.bam
    jobname="$short.CNV.FACETS"
    sbatchname="$short.CNV.FACETS.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=default1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=160G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    echo "mamba activate py27" >>$sbatchname
    echo "$cnvfacets -t $fileT -n $fileN -vcf $dbsnppath -o $short.facet_out --snp-nprocs $ncpu" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 6.1phyloWGS (top-level): clone evolution from SSMs + allele-specific CNVs
#
# Two environments: the FACETS->phyloWGS converter and the top-5 selector run
# under the base python3, while the phyloWGS tools need the py27 env (which
# ships no python3, so PY3 is captured before activation).
###############################################################################
gen_phylowgs_batcher() {
    cat <<'EOF'
#!/bin/bash
ncpu=8
for i in $(ls *.facet_out.vcf.gz); do
    short=$(echo $i | sed 's/\(.*\)\.facet_out\.vcf\.gz/\1/')
    jobname="$short.phyloWGS"
    sbatchname="$short.phyloWGS.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=default1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=64G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    echo "PY3=\$(command -v python3)" >>$sbatchname
    # 1. FACETS segments -> phyloWGS CNV table (base python3)
    echo "\"\$PY3\" facets_vcf_to_phylowgs_cnvs.py --facets $short.facet_out.vcf.gz --out $short.cnvs.txt" >>$sbatchname
    # 2. Build phyloWGS inputs (py27)
    echo "mamba activate py27" >>$sbatchname
    echo "create_phylowgs_inputs.py --vcf-type S1=mutect_smchet --tumor-sample PC_$short.Tumor --cnvs S1=$short.cnvs.txt --regions all --output-variants $short.ssm_data.txt --output-cnvs $short.cnv_data.txt --output-params $short.params.json S1=$short.mutect2.filtered.PASS.norm.snpeff.missense_frameshift.vcf" >>$sbatchname
    # 3. MCMC -- 4 chains run as parallel subprocesses
    echo "multievolve.py --num-chains 5 --ssms $short.ssm_data.txt --cnvs $short.cnv_data.txt --params $short.params.json --output-dir $short.smoke --burnin-samples 1000 --mcmc-samples 2500" >>$sbatchname
    # 4. Write trees
    echo "write_results.py --include-ssm-names $short $short.smoke/trees.zip $short.summ.json.gz $short.muts.json.gz $short.mutass.zip" >>$sbatchname
    # 5. Top-5 trees by log-likelihood (base python3, script carries its own)
    echo "bash selecttop5tree.sh $short" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 7.HLA_type (top-level): HLA-LA germline HLA typing on the blood BAM
###############################################################################
gen_hla_batcher() {
    cat <<'EOF'
#!/bin/bash
ncpu=16
for i in $(ls *.markdup.recal.bam); do
    short=$(echo $i | sed 's/\(.*\)\.Blood\.markdup\.recal\.bam/\1/')
    jobname="$short.HLA-LA"
    sbatchname="$short.HLA-LA.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=gpu1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=100G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    echo "mamba activate HLA-LA" >>$sbatchname
    echo "unset PERL5LIB" >>$sbatchname
    echo "HLA-LA.pl --workingDir ./ --BAM $i --sampleID $short --maxThreads $ncpu --graph PRG_MHC_GRCh38_withIMGT" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# 8.VacQ (top-level): VEP annotation + pVACseq neoantigen prediction
#
# VEP runs first, producing the .vep.vcf.gz that pVACseq consumes. The HLA
# allele list and the --tumor-purity argument are derived per sample by
# vacq_prepare.py under the base python3; the pvactools env is activated only
# around VEP and pVACseq themselves.
###############################################################################
gen_vacq_batcher() {
    cat <<'EOF'
#!/bin/bash
vepcache="/sibcb1/jizhelab1/hukaining/Biosoft/VEP/cache/"
vepplugins="/sibcb1/jizhelab1/hukaining/Biosoft/VEP/VEP_plugins/"
iedbdir="/sibcb1/jizhelab1/hukaining/Biosoft/forpVACtools/"
neodir="/sibcb1/jizhelab1/hukaining/202608/OE_WGS/NeoantigenEditing-main"
blastp_exe="/sibcb1/jizhelab1/hukaining/miniforge3/bin/blastp"
makeblastdb_exe="/sibcb1/jizhelab1/hukaining/miniforge3/bin/makeblastdb"
reffa="/sibcb1/jizhelab1/hukaining/genome/hg38/GRCh38.primary_assembly.genome.fa"
ncpu=16
for i in $(ls *.mutect2_strelka.both.PASS.snpeff.missense_frameshift.vcf); do
    stem=$(echo $i | sed 's/\.vcf$//')
    short=$(echo $stem | sed 's/\(.*\)\.mutect2_strelka\.both\.PASS\.snpeff\.missense_frameshift/\1/')
    pvout="pvacseq_output/MHC_Class_I/PC_$short.Tumor.MHC_I"
    jobname="$short.VacQ"
    sbatchname="$short.VacQ.sbatch"
    echo "#!/bin/bash" >$sbatchname
    echo "#SBATCH --job-name=$jobname" >>$sbatchname
    echo "#SBATCH --output=$jobname.out" >>$sbatchname
    echo "#SBATCH --error=$jobname.err" >>$sbatchname
    echo "#SBATCH --partition=default1" >>$sbatchname
    echo "#SBATCH --ntasks=1" >>$sbatchname
    echo "#SBATCH --cpus-per-task=$ncpu" >>$sbatchname
    echo "#SBATCH --mem=120G" >>$sbatchname
    echo "source ~/.bashrc" >>$sbatchname
    echo "PY3=\$(command -v python3)" >>$sbatchname
    # Confident class I HLA alleles and the purity flag (base python3)
    echo "HLA_ALLELES=\$(\"\$PY3\" vacq_prepare.py hla R1_bestguess_G.txt)" >>$sbatchname
    echo "PURITY_ARG=\$(\"\$PY3\" vacq_prepare.py purity $short.facet_out.vcf.gz)" >>$sbatchname
    # ---- Pick the VEP input ----
    # Prefer the both-callers-PASS chain, but fall back to Mutect2-only: the
    # snpEff protein_coding filter can remove every site of the both-PASS chain
    # (BEG091), and VEP fails outright on a header-only VCF.
    echo "VCF_BOTH=$short.mutect2_strelka.both.PASS.snpeff.missense_frameshift.vcf" >>$sbatchname
    echo "VCF_M2=$short.mutect2.filtered.PASS.norm.snpeff.missense_frameshift.vcf" >>$sbatchname
    echo "INPUT=''" >>$sbatchname
    echo "N_INPUT=0" >>$sbatchname
    echo "for cand in \"\$VCF_BOTH\" \"\$VCF_M2\"; do" >>$sbatchname
    echo "    if [ -f \"\$cand\" ]; then" >>$sbatchname
    echo "        n=\$(grep -vc '^#' \"\$cand\" || true)" >>$sbatchname
    echo "        if [ \"\$n\" -gt 0 ]; then INPUT=\"\$cand\"; N_INPUT=\$n; break; fi" >>$sbatchname
    echo "        echo \"warning: \$cand has 0 variants\" >&2" >>$sbatchname
    echo "    fi" >>$sbatchname
    echo "done" >>$sbatchname
    echo "if [ -z \"\$INPUT\" ]; then" >>$sbatchname
    echo "    echo \"ERROR: no missense/frameshift variants in either chain for $short; skipping VEP/pvacseq/quality scoring\" >&2" >>$sbatchname
    echo "    exit 0" >>$sbatchname
    echo "fi" >>$sbatchname
    echo "STEM=\${INPUT%.vcf}" >>$sbatchname
    echo "echo \"VEP input: \$INPUT (\$N_INPUT variants)\" >&2" >>$sbatchname
    # VEP annotation (pvactools env)
    echo "mamba activate pvactools" >>$sbatchname
    echo "vep -i \"\$STEM.vcf\" -o \"\$STEM.vep.vcf\" --vcf --assembly GRCh38 --cache --offline --dir_cache $vepcache --dir_plugins $vepplugins --fasta $reffa --symbol --terms SO --tsl --mane_select --canonical --biotype -hgvs --plugin Wildtype --plugin Frameshift --fork 8 --buffer_size 5000 --no_stats --force_overwrite" >>$sbatchname
    echo "bgzip -f \"\$STEM.vep.vcf\"" >>$sbatchname
    echo "tabix -f -p vcf \"\$STEM.vep.vcf.gz\"" >>$sbatchname
    # Neoantigen prediction
    echo "pvacseq run \"\$STEM.vep.vcf.gz\" PC_$short.Tumor \"\$HLA_ALLELES\" NetMHCpan pvacseq_output -e1 9 --normal-sample-name $short.Blood --iedb-install-directory $iedbdir --binding-threshold 4000 --n-threads $ncpu --tdna-vaf 0.02 --normal-vaf 0.04 --tdna-cov 5 \$PURITY_ARG" >>$sbatchname
    # The scoring steps consume pVACseq output; stop early rather than letting
    # both scripts fail on a missing file.
    echo "if [ ! -f $pvout.all_epitopes.tsv ]; then" >>$sbatchname
    echo "    echo \"ERROR: pvacseq produced no output ($pvout.all_epitopes.tsv); skipping quality scoring\" >&2" >>$sbatchname
    echo "    exit 0" >>$sbatchname
    echo "fi" >>$sbatchname
    # Map variants to phyloWGS top-tree clones and estimate CCF.
    # CCF is prevalence/purity, so it needs BOTH a usable purity and the
    # phyloWGS top trees. When either is unavailable the mapping is skipped and
    # the quality scoring below still runs, just without clonality columns.
    echo "PHYLOWGS_OK=1" >>$sbatchname
    echo "for req in $short.ssm_data.txt $short.top5.summ.json.gz $short.top5.mutass.zip $short.top5_tree_ids.tsv; do" >>$sbatchname
    echo "    if [ ! -f \"\$req\" ]; then echo \"missing phyloWGS input: \$req\" >&2; PHYLOWGS_OK=0; fi" >>$sbatchname
    echo "done" >>$sbatchname
    echo "if [ -n \"\$PURITY_ARG\" ] && [ \"\$PHYLOWGS_OK\" = 1 ]; then" >>$sbatchname
    echo "    PURITY=\"\${PURITY_ARG#--tumor-purity }\"" >>$sbatchname
    echo "    \"\$PY3\" phylowgs_to_clonality.py --ssm-data $short.ssm_data.txt --summary-gz $short.top5.summ.json.gz --mutass-zip $short.top5.mutass.zip --pvac-tsv $pvout.all_epitopes.tsv --tree-ids $short.top5_tree_ids.tsv --tumor-purity \"\$PURITY\" --output $short.ccf_clonality.tsv" >>$sbatchname
    echo "    CLONALITY_ARG=\"--clonality-tsv $short.ccf_clonality.tsv\"" >>$sbatchname
    echo "else" >>$sbatchname
    echo "    echo 'skipping CCF/clonality mapping (needs a usable purity AND phyloWGS top trees); quality scoring continues without it' >&2" >>$sbatchname
    echo "    CLONALITY_ARG=''" >>$sbatchname
    echo "fi" >>$sbatchname
    # Neoantigen quality: Q2022 (NeoantigenEditing) and Q2023 9-mer
    echo "\"\$PY3\" convert_pvacseq_quality.py --pvac-tsv $pvout.all_epitopes.tsv --filtered-tsv $pvout.filtered.tsv --iedb-fasta $neodir/data/iedb.fasta --model-json $neodir/data/epitope_distance_model_parameters.json --sample-id $short --output-dir quality_output --blastp $blastp_exe --makeblastdb $makeblastdb_exe \$CLONALITY_ARG" >>$sbatchname
    echo "sbatch $sbatchname"
done
EOF
}

###############################################################################
# Main
###############################################################################
echo "=== OE_WGS Downstream Pipeline Setup (stages 3.1-7) ==="
echo "Samples: ${ALL_SAMPLES[*]}"
echo ""

for sample in "${ALL_SAMPLES[@]}"; do
    SAMPLE_DIR="${BASE_DIR}/${sample}"

    if [[ ! -d "$SAMPLE_DIR" ]]; then
        echo "WARNING: ${SAMPLE_DIR} does not exist, skipping. Run setup_pipeline.sh first."
        continue
    fi

    echo "--- Setting up downstream stages for ${sample} ---"

    # ---- Create directories ----
    mkdir -p "${SAMPLE_DIR}/WGS/3.1Pileup_Conta"
    mkdir -p "${SAMPLE_DIR}/WGS_blood/3.1Pileup"
    mkdir -p "${SAMPLE_DIR}/3.2CalCont"
    mkdir -p "${SAMPLE_DIR}/4.SNV"
    mkdir -p "${SAMPLE_DIR}/4_1.St_SNV"
    mkdir -p "${SAMPLE_DIR}/5.FilterSNV"
    mkdir -p "${SAMPLE_DIR}/6.CNV"
    mkdir -p "${SAMPLE_DIR}/6.1phyloWGS"
    mkdir -p "${SAMPLE_DIR}/7.HLA_type"
    mkdir -p "${SAMPLE_DIR}/8.VacQ"

    # ---- Write batcher scripts ----
    gen_pileup_tumor_batcher > "${SAMPLE_DIR}/WGS/3.1Pileup_Conta/batPileup_Conta.sh"
    gen_pileup_blood_batcher > "${SAMPLE_DIR}/WGS_blood/3.1Pileup/batPileup.sh"
    gen_calcont_batcher     > "${SAMPLE_DIR}/3.2CalCont/batCalCont.sh"
    gen_mutect2_batcher     > "${SAMPLE_DIR}/4.SNV/batMutect2_BQSR.sh"
    gen_f1r2_batcher        > "${SAMPLE_DIR}/4.SNV/batF1R2Model.sh"
    gen_strelka_batcher     > "${SAMPLE_DIR}/4_1.St_SNV/batStelka_BQSR.sh"
    gen_filter_batcher      > "${SAMPLE_DIR}/5.FilterSNV/batFilter.sh"
    gen_cnv_batcher         > "${SAMPLE_DIR}/6.CNV/batCNV_FACETS.sh"
    gen_phylowgs_batcher    > "${SAMPLE_DIR}/6.1phyloWGS/batPhyloWGS.sh"
    gen_hla_batcher         > "${SAMPLE_DIR}/7.HLA_type/batHLAtp.sh"
    gen_vacq_batcher        > "${SAMPLE_DIR}/8.VacQ/batVacQ.sh"

    # ---- Copy the shared phyloWGS helpers (masters live in 6.1phyloWGS/) ----
    cp -p "${BASE_DIR}/6.1phyloWGS/facets_vcf_to_phylowgs_cnvs.py" "${SAMPLE_DIR}/6.1phyloWGS/"
    cp -p "${BASE_DIR}/6.1phyloWGS/selecttop5tree.sh" "${SAMPLE_DIR}/6.1phyloWGS/"

    # ---- Copy the shared VacQ helpers (masters live in 8.VacQ/) ----
    cp -p "${BASE_DIR}/8.VacQ/vacq_prepare.py" "${SAMPLE_DIR}/8.VacQ/"
    cp -p "${BASE_DIR}/8.VacQ/phylowgs_to_clonality.py" "${SAMPLE_DIR}/8.VacQ/"
    cp -p "${BASE_DIR}/8.VacQ/convert_pvacseq_quality.py" "${SAMPLE_DIR}/8.VacQ/"

    chmod +x "${SAMPLE_DIR}/WGS/3.1Pileup_Conta/batPileup_Conta.sh"
    chmod +x "${SAMPLE_DIR}/WGS_blood/3.1Pileup/batPileup.sh"
    chmod +x "${SAMPLE_DIR}/3.2CalCont/batCalCont.sh"
    chmod +x "${SAMPLE_DIR}/4.SNV/batMutect2_BQSR.sh"
    chmod +x "${SAMPLE_DIR}/4.SNV/batF1R2Model.sh"
    chmod +x "${SAMPLE_DIR}/4_1.St_SNV/batStelka_BQSR.sh"
    chmod +x "${SAMPLE_DIR}/5.FilterSNV/batFilter.sh"
    chmod +x "${SAMPLE_DIR}/6.CNV/batCNV_FACETS.sh"
    chmod +x "${SAMPLE_DIR}/6.1phyloWGS/batPhyloWGS.sh"
    chmod +x "${SAMPLE_DIR}/6.1phyloWGS/selecttop5tree.sh"
    chmod +x "${SAMPLE_DIR}/7.HLA_type/batHLAtp.sh"
    chmod +x "${SAMPLE_DIR}/8.VacQ/batVacQ.sh"

    echo "  Created: ${SAMPLE_DIR}/WGS/3.1Pileup_Conta"
    echo "  Created: ${SAMPLE_DIR}/WGS_blood/3.1Pileup"
    echo "  Created: ${SAMPLE_DIR}/3.2CalCont"
    echo "  Created: ${SAMPLE_DIR}/4.SNV (Mutect2 + F1R2 model)"
    echo "  Created: ${SAMPLE_DIR}/4_1.St_SNV"
    echo "  Created: ${SAMPLE_DIR}/5.FilterSNV"
    echo "  Created: ${SAMPLE_DIR}/6.CNV (FACETS)"
    echo "  Created: ${SAMPLE_DIR}/6.1phyloWGS (+ helper scripts copied)"
    echo "  Created: ${SAMPLE_DIR}/7.HLA_type"
    echo "  Created: ${SAMPLE_DIR}/8.VacQ (+ helper script copied)"
    echo ""
done

echo "=== Done ==="
echo "Downstream stages created for ${#ALL_SAMPLES[@]} samples."
echo ""

# Generate the pipeline dependency map for reference
cat <<'MAP_EOF'

Dependency graph (after steps 1-3 complete, in order):

  WGS/3.BQSR ──────────────────┬──► WGS/3.1Pileup_Conta ──► 3.2CalCont* ──┐
  (tumor recal.bam)            │                                            │
                               ├──► 4.SNV* ──► 4.SNV/F1R2* ─────────────────┼──► 5.FilterSNV*
                               ├──► 4_1.St_SNV*                            │
                               │                                            │
  WGS_blood/3.BQSR ────────────┼──► WGS_blood/3.1Pileup ────► 3.2CalCont*──┘
  (blood recal.bam)            │
                               ├──► 4.SNV*
                               ├──► 4_1.St_SNV*
                               ├──► 6.CNV*         (FACETS, both BAMs)
                               └──► 7.HLA_type*   (blood BAM only)

  5.FilterSNV ─┐
               ├──► 6.1phyloWGS*   (needs both the SSM VCF and 6.CNV output)
  6.CNV ───────┘

  5.FilterSNV (both-PASS missense/frameshift VCF) ─┐
  6.CNV (facet_out.vcf.gz, for purity)            ├──► 8.VacQ*
  7.HLA_type (R1_bestguess_G.txt) ────────────────┘

  (* = top-level stage, reads from both WGS/ and WGS_blood/)

Run order:
  1. WGS/3.1Pileup_Conta + WGS_blood/3.1Pileup  (parallel)
     bash link_stage_outputs.sh <sample> 3.1     (wire pileup tables)
  2. 3.2CalCont
     bash link_stage_outputs.sh <sample> 3.2     (wire contamination tables)
  3. 4.SNV/Mutect2 + 4_1.St_SNV + 6.CNV + 7.HLA_type  (parallel)
     bash link_stage_outputs.sh <sample> 4
  4. 4.SNV/F1R2 model (LearnReadOrientationModel, depends on Mutect2 f1r2.tar.gz)
     bash link_stage_outputs.sh <sample> 4.f1r2   (wire read-orientation-model)
  5. 5.FilterSNV
     (reads VCF + stats from step 3, contamination from step 2,
      and orientation model from step 4)
  6. 6.1phyloWGS (AFTER both 5.FilterSNV and 6.CNV finish)
     bash link_stage_outputs.sh <sample> 6.1     (wire SSM VCF + FACETS VCF)
  7. 8.VacQ (AFTER 5.FilterSNV, 6.CNV and 7.HLA_type finish)
     bash link_stage_outputs.sh <sample> 8       (wire VCF + FACETS VCF + HLA)
MAP_EOF
