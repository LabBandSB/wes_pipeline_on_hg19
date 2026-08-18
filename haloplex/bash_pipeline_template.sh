# HaloPlex Cardio panel pipeline — step-token template
#
# Derived from ../new_version/bash_pipeline_template.sh (WES).
# Every deviation from the WES template is marked with a "HALOPLEX:" comment
# explaining WHY. Read haloplex/README.md before changing any of them.

mkdir -p ${alignment_dir}
# HALOPLEX: the panel is 51.9 kb, not a 71 Mb exome. 64G of heap was sized for
# whole-exome BAMs; here it only slows JVM startup. 8G is generous for this data.
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
[ -f ${read1} ] && \
[ -f ${read2} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${bwa} mem -M -t ${threads} ${ref} ${read1} ${read2} > ${output_file} && \
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
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${samtools} view -bT ${ref} ${input_file} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# samtools sort bam to bam
token="${alignment_dir}/token.${sample}.bam_2_bam_samtools_sort"
input_file="${alignment_dir}/${sample}.samtools_view.bam"
output_file="${alignment_dir}/${sample}.samtools_sort.bam"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${samtools} sort -l 9 -O bam -T ${alignment_dir}/${sample}.sorted.tmp ${input_file} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# picard_ARRG
token="${alignment_dir}/token.${sample}.bam_2_bam_picard_ARRG"
input_file="${alignment_dir}/${sample}.samtools_sort.bam"
output_file="${alignment_dir}/${sample}.picard_ARRG.bam"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
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
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


###################################################################
# HALOPLEX: MarkDuplicates is DELIBERATELY ABSENT.
#
# The WES template runs picard MarkDuplicates here. For HaloPlex that step
# does not remove PCR duplicates — it removes real, independent molecules.
#
# HaloPlex fragments DNA with restriction enzymes, not by random shearing, so
# every read covering a given amplicon starts and ends at the SAME coordinates
# by construction. Positional deduplication assumes coordinate collisions are
# improbable and therefore indicate PCR copies. That assumption is false here.
# Running MarkDuplicates on this data typically flags 70-90% of reads and
# collapses coverage to near-nothing.
#
# Rule of thumb that generalises: positional dedup is valid only when
# fragmentation is RANDOM (hybrid capture, WGS). It is invalid for every
# amplicon chemistry — HaloPlex, AmpliSeq, QIAseq, Fluidigm.
#
# The exception is HaloPlexHS, which carries UMIs. That data MUST be
# deduplicated, but by molecular barcode (Agilent AGeNT LocatIt), never by
# coordinate. These samples are classic HaloPlex — no UMIs — so nothing to do.
###################################################################


###################################################################
# HALOPLEX: BQSR (BaseRecalibrator / AnalyzeCovariates / PrintReads) is
# DELIBERATELY ABSENT.
#
# The WES template restricts BaseRecalibrator with -L ${target_region}. With a
# 71 Mb exome that leaves plenty of known sites to build covariate tables from.
# This panel is 51.9 kb across 60 regions — roughly 1400x less territory. The
# recalibration model would be fitted on a few hundred known sites, which is far
# below what GATK needs for a stable estimate; the result is either an error or,
# worse, a confidently wrong recalibration.
#
# Running BQSR unrestricted does not help either: panel data has very few
# off-target reads to widen the training set with.
#
# Skipping BQSR on small panels is the standard, documented choice. Base
# qualities from the sequencer are used as-is.
###################################################################


# #######
# HALOPLEX: HaplotypeCaller straight off the aligned BAM (no MD, no BQSR).
# -L restricts calling to the panel; -ip 50 pads each region by 50 bp because
# HaloPlex amplicon ends are fixed and variants sitting right on a boundary
# would otherwise be missed or called with truncated context.
token="${alignment_dir}/token.${sample}.bam_2_vcf_gatk_HC"
input_file="${alignment_dir}/${sample}.picard_ARRG.bam"
output_file="${alignment_dir}/${sample}.gatk_HC.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -Xmx${XMXVALUE} -T HaplotypeCaller \
  -R ${ref} \
  -I ${input_file} \
  --dbsnp ${dbsnp} \
  -L ${target_region} \
  -ip ${interval_padding} \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30 \
  -o ${output_file} \
  -nct ${threads} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


###################################################################
# HALOPLEX: VQSR (VariantRecalibrator / ApplyRecalibration) is DELIBERATELY
# ABSENT — replaced by hard filtering below.
#
# GATK's own requirement is at least 30 exomes or 1 whole genome to train
# VariantRecalibrator. The WES template already acknowledges this in a comment
# and works around it with --maxGaussians 1, which suppresses the "no data"
# crash without making the model meaningful.
#
# One sample over 51.9 kb yields on the order of tens of variants. There is no
# amount of tuning that makes a Gaussian mixture model trainable on that.
# GATK's documented fallback for small targets is hard filtering, which is what
# the SNP/INDEL VariantFiltration steps below do — and which the WES template
# already contained downstream of VQSR anyway.
#
# If we later switch to joint genotyping across all 143 samples, VQSR becomes
# worth revisiting — see README.md, "Open questions".
###################################################################


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_SV_SNP"
input_file="${alignment_dir}/${sample}.gatk_HC.vcf"
output_file="${alignment_dir}/${sample}.gatk_SV_SNP.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -Xmx${XMXVALUE} -T SelectVariants \
  -R ${ref} \
  -V ${input_file} \
  -selectType SNP \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
# HALOPLEX: MQRankSum and ReadPosRankSum are kept, but note they are computed
# from read-position distributions. With fixed amplicon ends these are less
# informative than on randomly sheared data — do not tighten them further
# without looking at the actual distributions first.
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_VF_SNP"
input_file="${alignment_dir}/${sample}.gatk_SV_SNP.vcf"
output_file="${alignment_dir}/${sample}.gatk_VF_SNP.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -Xmx${XMXVALUE} -T VariantFiltration \
  -R ${ref} \
  -V ${input_file} \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filterName "SNP_FAIL" \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_SV_INDEL"
input_file="${alignment_dir}/${sample}.gatk_HC.vcf"
output_file="${alignment_dir}/${sample}.gatk_SV_INDEL.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -Xmx${XMXVALUE} -T SelectVariants \
  -R ${ref} \
  -V ${input_file} \
  -selectType INDEL \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_VF_INDEL"
input_file="${alignment_dir}/${sample}.gatk_SV_INDEL.vcf"
output_file="${alignment_dir}/${sample}.gatk_VF_INDEL.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -Xmx${XMXVALUE} -T VariantFiltration \
  -R ${ref} \
  -V ${input_file} \
  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
  --filterName "INDEL_FAIL" \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_vcf_2_vcf_vcf_concat"
input_file="${alignment_dir}/${sample}.gatk_VF_SNP.vcf"
input_file_2="${alignment_dir}/${sample}.gatk_VF_INDEL.vcf"
output_file="${alignment_dir}/${sample}.vcf_concat.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
[ -f ${input_file_2} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${vcf_concat} ${input_file} ${input_file_2} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_vcf_sort"
input_file="${alignment_dir}/${sample}.vcf_concat.vcf"
output_file="${alignment_dir}/${sample}.vcf_sort.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
cat ${input_file} | ${vcf_sort} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_FINAL"
input_file="${alignment_dir}/${sample}.vcf_sort.vcf"
output_file="${alignment_dir}/${sample}.FINAL.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
cp ${input_file} ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


###################################################################
# HALOPLEX: panel coverage QC — NEW step, absent from the WES template.
#
# On a 51.9 kb panel, coverage is the first thing that goes wrong and the
# easiest to miss. A sample that failed enrichment still produces a valid-looking
# VCF, just with fewer calls. This step makes that visible per sample.
#
# Note the two different BEDs: depth is measured over the Amplicons BED (1202
# amplicons, 240 kb) because that is what was physically amplified, while
# calling above used the Covered BED (60 regions, 51.9 kb). Using the wrong one
# is the classic HaloPlex QC mistake.
###################################################################
token="${alignment_dir}/token.${sample}.bam_2_txt_coverage_qc"
input_file="${alignment_dir}/${sample}.picard_ARRG.bam"
output_file="${alignment_dir}/${sample}.panel_coverage.txt"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
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
# HALOPLEX: ANNOVAR annotation — NEW step, absent from the WES template.
# Databases live on Kilo, see README.md. Uses hg19 to match the panel design.
###################################################################
token="${alignment_dir}/token.${sample}.vcf_2_annovar"
input_file="${alignment_dir}/${sample}.FINAL.vcf"
output_file="${alignment_dir}/${sample}.annovar.hg19_multianno.txt"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
[ -d ${annovar_humandb} ] && \
rm -f ${output_file} && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
perl ${annovar_dir}/table_annovar.pl ${input_file} ${annovar_humandb} \
  --buildver hg19 \
  --outfile ${alignment_dir}/${sample}.annovar \
  --remove \
  --protocol ${annovar_protocol} \
  --operation ${annovar_operation} \
  --nastring . \
  --vcfinput && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
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


echo "!!! HALOPLEX PIPELINE DONE FOR SAMPLE=${sample} !!!"
