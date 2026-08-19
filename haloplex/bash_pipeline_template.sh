# HaloPlex Cardio panel — PER-SAMPLE stage, step-token template.
#
# Produces one GVCF per sample. Genotyping happens later, once, across the whole
# cohort — see run_cohort.sh. Nothing here emits a final VCF on its own.
#
# Derived from ../new_version/bash_pipeline_template.sh (WES). Every deviation is
# marked "HALOPLEX:" with the reason. Read haloplex/README.md before changing any.

###################################################################
# HALOPLEX: every step tests its input with -s (non-empty), not -f (exists),
# and verifies its own output is non-empty before writing the success token.
#
# Why this is not paranoia. On 2026-08-19 a `bwa mem` process was killed while
# its wrapper script kept running. bwa left a 0-byte SAM. `samtools view` on an
# empty SAM exits 0 and writes an empty BAM — so the token was written. Sort,
# AddOrReplaceReadGroups and the coverage step all then "succeeded" on empty
# input, the FINAL_LOCK was set, and Haloplex_472 looked complete in six seconds
# with 0 of its 2 432 860 read pairs. Nothing errored. A re-run would have
# skipped it.
#
# Emptiness must not propagate as success.
###################################################################

mkdir -p ${alignment_dir}
# HALOPLEX: the panel is 464 kb, not a 71 Mb exome. 64G of heap was sized for
# whole-exome BAMs; here it only slows JVM startup.
XMXVALUE="8G"


###################################################################
# first and final locks declaration BEGIN
###################################################################
FIRST_LOCK="${alignment_dir}/token.${sample}.__FIRST_LOCK__"
FINAL_LOCK="${alignment_dir}/token.${sample}.__FINAL_LOCK__"
HISTORY_LOCK="${alignment_dir}/token.${sample}._HISTORY_LOCK_"

# tokens mode description:
# (x, y) - first and final locks mode
# (0, 0) - not started yet - feel free to start
# (1, 0) - started, but not ended - don't touch - it runs somewhere
# (0, 1) - don't started, but already ended - how? why ? - for now treated as (0, 0)
# (1, 1) - started and ended - you can try rerun

# (1, 1) - check if already completed to rerun
[ -f  ${FINAL_LOCK} ] && \
rm -f ${FIRST_LOCK} && \
rm -f ${FINAL_LOCK} && \
echo "RERUN ${sample}"

# (0, 0), (0, 1) - lock sample to prevent runnig from 2 servers
[ ! -f ${FIRST_LOCK} ] && \
dt1dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1dt1} > ${FIRST_LOCK} && \
echo ${dt1dt1} >> ${HISTORY_LOCK} && \
rm -f ${FINAL_LOCK} \
|| exit 1
###################################################################
# first and final locks declaration END
###################################################################


# bwa alignment
token="${alignment_dir}/token.${sample}.fastq_2_sam_bwa_mem"
output_file="${alignment_dir}/${sample}.bwa_mem.sam"
[ ! -f ${token} ] && \
[ -s ${read1} ] && \
[ -s ${read2} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${bwa} mem -M -t ${threads} ${ref} ${read1} ${read2} > ${output_file} && \
[ -s ${output_file} ] && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# samtools convert sam to bam
token="${alignment_dir}/token.${sample}.sam_2_bam_samtools_view"
input_file="${alignment_dir}/${sample}.bwa_mem.sam"
output_file="${alignment_dir}/${sample}.samtools_view.bam"
[ ! -f ${token} ] && \
[ -s ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${samtools} view -bT ${ref} ${input_file} > ${output_file} && \
[ -s ${output_file} ] && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# samtools sort bam to bam
# HALOPLEX: -l 6 (default) instead of -l 9, and -@ threads instead of single-core.
# Measured on this cohort: sort at -l 9 single-threaded had a median of 904 s and a
# worst case of 97 min — 58% of total per-sample time, more than bwa and
# HaplotypeCaller combined. Level 9 buys roughly a tenth of the file size for
# several times the CPU, and the sort was using one core while fifteen sat idle.
token="${alignment_dir}/token.${sample}.bam_2_bam_samtools_sort"
input_file="${alignment_dir}/${sample}.samtools_view.bam"
output_file="${alignment_dir}/${sample}.samtools_sort.bam"
[ ! -f ${token} ] && \
[ -s ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${samtools} sort -l 6 -@ ${threads} -O bam -T ${alignment_dir}/${sample}.sorted.tmp ${input_file} > ${output_file} && \
[ -s ${output_file} ] && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# picard AddOrReplaceReadGroups
# HALOPLEX: RGSM is what GenotypeGVCFs uses as the column name in the cohort VCF.
# If it is wrong here, it is wrong in the final matrix and nowhere else will tell you.
token="${alignment_dir}/token.${sample}.bam_2_bam_picard_ARRG"
input_file="${alignment_dir}/${sample}.samtools_sort.bam"
output_file="${alignment_dir}/${sample}.picard_ARRG.bam"
[ ! -f ${token} ] && \
[ -s ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${picard} AddOrReplaceReadGroups \
  INPUT=${input_file} \
  OUTPUT=${output_file} \
  SORT_ORDER=coordinate \
  RGID=${sample} \
  RGLB=${sample} \
  RGPL=ILLUMINA \
  RGPU=${platform_unit} \
  RGSM=${sample} \
  RGCN=NLA \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_RECORDS_IN_RAM=1000000 && \
[ "$(${samtools} view -c ${output_file})" -gt 0 ] && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


###################################################################
# HALOPLEX: MarkDuplicates is DELIBERATELY ABSENT.
#
# The WES template runs picard MarkDuplicates here. For HaloPlex that step does
# not remove PCR duplicates — it removes real, independent molecules.
#
# HaloPlex fragments DNA with restriction enzymes, not by random shearing, so
# every read covering a given amplicon starts and ends at the SAME coordinates
# by construction. Positional deduplication assumes coordinate collisions are
# improbable and therefore indicate PCR copies. That assumption is false here.
# Running MarkDuplicates on this data typically flags 70-90% of reads and
# collapses coverage to near-nothing.
#
# Rule that generalises: positional dedup is valid only when fragmentation is
# RANDOM (hybrid capture, WGS). It is invalid for every amplicon chemistry —
# HaloPlex, AmpliSeq, QIAseq, Fluidigm.
#
# The exception is HaloPlexHS, which carries UMIs. That data MUST be
# deduplicated, but by molecular barcode (Agilent AGeNT LocatIt), never by
# coordinate. These samples are classic HaloPlex — no UMIs — so nothing to do.
###################################################################


###################################################################
# HALOPLEX: BQSR is DELIBERATELY ABSENT.
#
# The WES template restricts BaseRecalibrator with -L ${target_region}. With a
# 71 Mb exome that leaves plenty of known sites to build covariate tables from.
# This panel is 464 kb across 2007 regions — roughly 150x less territory. The
# model would be fitted on a few thousand known sites at best, below what GATK
# needs for a stable estimate; the result is either an error or, worse, a
# confidently wrong recalibration applied to every base.
#
# Running BQSR unrestricted does not help: panel data has very few off-target
# reads to widen the training set with. Base qualities are used as-is.
###################################################################


# #######
# HALOPLEX: HaplotypeCaller in GVCF mode, straight off the aligned BAM.
#
#   -ERC GVCF  emits a per-sample GVCF with reference blocks instead of a VCF.
#              Genotyping is deferred to the cohort stage (run_cohort.sh), which
#              is what makes a squared-off matrix possible: every sample gets a
#              genotype at every site the cohort varies at, and "homozygous
#              reference" becomes distinguishable from "no coverage".
#   -L / -ip   restrict to the panel, padded by 50 bp — amplicon ends are fixed,
#              so a variant on a region boundary would otherwise be called with
#              truncated context or missed.
token="${alignment_dir}/token.${sample}.bam_2_gvcf_gatk_HC"
input_file="${alignment_dir}/${sample}.picard_ARRG.bam"
output_file="${alignment_dir}/${sample}.g.vcf.gz"
[ ! -f ${token} ] && \
[ -s ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} --java-options "-Xmx${XMXVALUE}" HaplotypeCaller \
  -R ${ref} \
  -I ${input_file} \
  -O ${output_file} \
  -L ${target_region} \
  -ip ${interval_padding} \
  --dbsnp ${dbsnp} \
  -ERC GVCF \
  --native-pair-hmm-threads ${threads} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


###################################################################
# HALOPLEX: panel coverage QC — NEW step, absent from the WES template.
#
# On a 464 kb panel, coverage is the first thing that goes wrong and the
# easiest to miss: a sample that failed enrichment still produces a valid-looking
# GVCF, just with fewer calls. This makes that visible per sample.
#
# Note the two different BEDs. Depth is measured over the Amplicons BED (1202
# amplicons, 240 kb) because that is what was physically amplified, while calling
# above used the Covered BED (60 regions, 464 kb). Using the wrong one is the
# classic HaloPlex QC mistake — it yields plausible numbers that mean something
# else entirely.
###################################################################
token="${alignment_dir}/token.${sample}.bam_2_txt_coverage_qc"
input_file="${alignment_dir}/${sample}.picard_ARRG.bam"
output_file="${alignment_dir}/${sample}.panel_coverage.txt"
[ ! -f ${token} ] && \
[ -s ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${samtools} depth -a -b ${amplicons_region} ${input_file} \
  | awk '{s+=$3; n++; if($3>=20) c20++; if($3>=100) c100++} \
         END {printf "sample\tmean_depth\tbases\tpct_ge20x\tpct_ge100x\n"; \
              printf "%s\t%.1f\t%d\t%.2f\t%.2f\n", "'"${sample}"'", s/n, n, 100*c20/n, 100*c100/n}' \
  > ${output_file} && \
cat ${output_file} && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


###################################################################
# final locks declaration BEGIN
###################################################################
dt2dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1dt1} ${dt2dt2} > ${FINAL_LOCK}
###################################################################
# final locks declaration END
###################################################################


echo "!!! PER-SAMPLE STAGE DONE FOR SAMPLE=${sample} -> ${alignment_dir}/${sample}.g.vcf.gz !!!"
echo "!!! run run_cohort.sh once every sample has finished !!!"
