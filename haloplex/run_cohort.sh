#!/bin/bash
# HaloPlex Cardio panel — COHORT stage.
#
# Run once, after every per-sample script has produced its <sample>.g.vcf.gz.
#
#   conda activate haloplex
#   ./run_cohort.sh settings.json
#
# Joint genotyping across all samples is what the per-sample GVCFs were for. It
# buys three things a per-sample pipeline cannot give:
#
#   1. A squared-off matrix. Every sample gets a genotype at every site the
#      cohort varies at, so "homozygous reference" stops being indistinguishable
#      from "no coverage here". For a cohort VCF that distinction is the point.
#   2. Sensitivity at low coverage. A weak signal in one sample is evaluated in
#      the light of confident calls at the same site in the other 142.
#   3. One consistent set of thresholds applied across the cohort, instead of
#      143 independent decisions.
#
# What it does NOT buy, despite a common expectation: VQSR. 143 samples over
# 51.9 kb yield on the order of hundreds to ~1500 variant sites. VariantRecalibrator
# needs thousands to fit a Gaussian mixture. Hard filtering below stays.

set -u
SETTINGS="${1:-settings.json}"
[ -f "$SETTINGS" ] || { echo "usage: $0 settings.json" >&2; exit 1; }

get() { python3 -c "import json,sys; print(json.load(open('$SETTINGS')).get('$1',''))"; }

REF=$(get ref)
DBSNP=$(get dbsnp)
TARGET=$(get target_region)
PADDING=$(get interval_padding)
PROJECT=$(get project_dir)
COHORT=$(get cohort_dir)
COHORT_NAME=$(get cohort_name)
THREADS=$(get threads)
ANNOVAR_DIR=$(get annovar_dir)
ANNOVAR_DB=$(get annovar_humandb)
ANNOVAR_PROTO=$(get annovar_protocol)
ANNOVAR_OPER=$(get annovar_operation)
XMX="${COHORT_XMX:-16G}"

mkdir -p "$COHORT"
LOG="$COHORT/_cohort.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== cohort stage started $(date) ==="

# ---------------------------------------------------------------- sample map
# GenomicsDBImport wants <sample-name>\t<path>. The sample name must match RGSM
# in the BAM, which the per-sample stage set to the directory name.
MAP="$COHORT/sample_map.tsv"
: > "$MAP"
n=0
for g in "$PROJECT"/*/*.g.vcf.gz; do
    [ -f "$g" ] || continue
    s=$(basename "$g" .g.vcf.gz)
    printf '%s\t%s\n' "$s" "$g" >> "$MAP"
    n=$((n+1))
done
echo "--- $n GVCFs found -> $MAP"
[ "$n" -gt 0 ] || { echo "ERROR: no GVCFs under $PROJECT/*/ — has the per-sample stage run?" >&2; exit 1; }

# Refuse to genotype a partial cohort silently. Joint calling on a subset gives
# different allele frequencies and different genotypes than the full set, so a
# half-finished run must not look like a finished one.
EXPECTED=$(ls -1 "$PROJECT" 2>/dev/null | wc -l)
if [ "$n" -lt "$EXPECTED" ]; then
    echo "WARNING: $n GVCFs but $EXPECTED sample directories — $((EXPECTED-n)) sample(s) unfinished."
    echo "         Joint genotyping a partial cohort changes allele frequencies and genotypes."
    echo "         Set COHORT_ALLOW_PARTIAL=1 to proceed anyway."
    [ "${COHORT_ALLOW_PARTIAL:-0}" = "1" ] || exit 1
fi

# ---------------------------------------------------------------- GenomicsDB
DB="$COHORT/genomicsdb"
if [ ! -d "$DB" ]; then
    echo "--- GenomicsDBImport $(date +%H:%M:%S)"
    gatk --java-options "-Xmx${XMX}" GenomicsDBImport \
        --sample-name-map "$MAP" \
        --genomicsdb-workspace-path "$DB" \
        -L "$TARGET" \
        --interval-padding "$PADDING" \
        --batch-size 50 \
        --reader-threads "$THREADS" || exit 1
else
    echo "--- GenomicsDBImport skipped (workspace exists: $DB)"
fi

# ---------------------------------------------------------------- genotyping
RAW="$COHORT/${COHORT_NAME}.raw.vcf.gz"
if [ ! -f "$RAW" ]; then
    echo "--- GenotypeGVCFs $(date +%H:%M:%S)"
    gatk --java-options "-Xmx${XMX}" GenotypeGVCFs \
        -R "$REF" \
        -V "gendb://$DB" \
        -O "$RAW" \
        --dbsnp "$DBSNP" \
        -L "$TARGET" \
        --interval-padding "$PADDING" || exit 1
else
    echo "--- GenotypeGVCFs skipped ($RAW exists)"
fi

# ---------------------------------------------------------------- hard filters
# GATK's documented hard-filter recommendations, applied separately to SNPs and
# indels because the sensible thresholds differ. See the note at the top of this
# file for why VQSR is not used instead.
#
# MQRankSum and ReadPosRankSum are computed from read-position distributions.
# With fixed amplicon ends these carry less information than on randomly sheared
# data — do not tighten them without first looking at the actual distributions.
SNP="$COHORT/${COHORT_NAME}.snp.filtered.vcf.gz"
if [ ! -f "$SNP" ]; then
    echo "--- SNPs: select + filter $(date +%H:%M:%S)"
    gatk --java-options "-Xmx${XMX}" SelectVariants \
        -R "$REF" -V "$RAW" -select-type SNP \
        -O "$COHORT/${COHORT_NAME}.snp.raw.vcf.gz" || exit 1
    gatk --java-options "-Xmx${XMX}" VariantFiltration \
        -R "$REF" -V "$COHORT/${COHORT_NAME}.snp.raw.vcf.gz" -O "$SNP" \
        --filter-expression "QD < 2.0"                --filter-name "QD2" \
        --filter-expression "FS > 60.0"               --filter-name "FS60" \
        --filter-expression "MQ < 40.0"               --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5"       --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0"   --filter-name "ReadPosRankSum-8" \
        --filter-expression "SOR > 3.0"               --filter-name "SOR3" || exit 1
fi

INDEL="$COHORT/${COHORT_NAME}.indel.filtered.vcf.gz"
if [ ! -f "$INDEL" ]; then
    echo "--- indels: select + filter $(date +%H:%M:%S)"
    gatk --java-options "-Xmx${XMX}" SelectVariants \
        -R "$REF" -V "$RAW" -select-type INDEL -select-type MIXED \
        -O "$COHORT/${COHORT_NAME}.indel.raw.vcf.gz" || exit 1
    gatk --java-options "-Xmx${XMX}" VariantFiltration \
        -R "$REF" -V "$COHORT/${COHORT_NAME}.indel.raw.vcf.gz" -O "$INDEL" \
        --filter-expression "QD < 2.0"                --filter-name "QD2" \
        --filter-expression "FS > 200.0"              --filter-name "FS200" \
        --filter-expression "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20" \
        --filter-expression "SOR > 10.0"              --filter-name "SOR10" || exit 1
fi

FINAL="$COHORT/${COHORT_NAME}.FINAL.vcf.gz"
if [ ! -f "$FINAL" ]; then
    echo "--- merge SNP + indel $(date +%H:%M:%S)"
    gatk --java-options "-Xmx${XMX}" MergeVcfs \
        -I "$SNP" -I "$INDEL" -O "$FINAL" || exit 1
fi

echo "--- cohort VCF: $FINAL"
bcftools index -f -t "$FINAL" 2>/dev/null
echo "    samples:  $(bcftools query -l "$FINAL" 2>/dev/null | wc -l)"
echo "    variants: $(bcftools index -n "$FINAL" 2>/dev/null)"
echo "    PASS:     $(bcftools view -f PASS -H "$FINAL" 2>/dev/null | wc -l)"

# ---------------------------------------------------------------- ANNOVAR
# Annotate once on the cohort VCF rather than 143 times per sample — same result,
# a fraction of the work, and one file to reason about.
if [ -d "$ANNOVAR_DB" ] && [ -f "$ANNOVAR_DIR/table_annovar.pl" ]; then
    OUT="$COHORT/${COHORT_NAME}.annovar"
    if [ ! -f "${OUT}.hg19_multianno.txt" ]; then
        echo "--- ANNOVAR $(date +%H:%M:%S)"
        perl "$ANNOVAR_DIR/table_annovar.pl" "$FINAL" "$ANNOVAR_DB" \
            --buildver hg19 \
            --outfile "$OUT" \
            --remove \
            --protocol "$ANNOVAR_PROTO" \
            --operation "$ANNOVAR_OPER" \
            --nastring . \
            --vcfinput || echo "ANNOVAR failed — cohort VCF is still valid"
    fi
    echo "    annotated: ${OUT}.hg19_multianno.txt"
else
    echo "--- ANNOVAR skipped: humandb or table_annovar.pl not reachable"
    echo "    ANNOVAR_DIR=$ANNOVAR_DIR"
    echo "    ANNOVAR_DB=$ANNOVAR_DB"
    echo "    (the NAS mounts are systemd automounts — 'ls' the path once to wake them)"
fi

# ---------------------------------------------------------------- coverage QC
QC="$COHORT/${COHORT_NAME}.panel_coverage.tsv"
echo -e "sample\tmean_depth\tbases\tpct_ge20x\tpct_ge100x" > "$QC"
for f in "$PROJECT"/*/*.panel_coverage.txt; do
    [ -f "$f" ] && tail -n +2 "$f" >> "$QC"
done
echo "--- coverage summary: $QC ($(( $(wc -l < "$QC") - 1 )) samples)"
awk -F'\t' 'NR>1 && $4+0 < 90 {print "    LOW COVERAGE: " $1 "  mean=" $2 "  %>=20x=" $4}' "$QC"

echo "=== cohort stage finished $(date) ==="
