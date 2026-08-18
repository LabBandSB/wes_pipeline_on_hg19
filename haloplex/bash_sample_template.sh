# sample aliases
read1="{read1}"
read2="{read2}"
sample="{sample}"
alignment_dir="{alignment_dir}"

# databases list
ref="{ref}"
dbsnp="{dbsnp}"
gold_indel="{gold_indel}"
hapmap_snp="{hapmap_snp}"
oneKG_indel="{oneKG_indel}"
oneKG_snp="{oneKG_snp}"
onmi_snp="{onmi_snp}"

# HALOPLEX: two different BEDs, used for two different things.
#   target_region    -> Covered.bed   (60 regions, 51.9 kb)  -> variant calling (-L)
#   amplicons_region -> Amplicons.bed (1202 amplicons, 240 kb) -> coverage QC
# Swapping them silently produces wrong results, not errors.
target_region="{target_region}"
amplicons_region="{amplicons_region}"
interval_padding="{interval_padding}"

# HALOPLEX: read-group platform unit. The WES template hardcoded SureSelectV4,
# which is a different capture chemistry — wrong metadata in every BAM header.
platform_unit="{platform_unit}"
threads="{threads}"

# annovar
annovar_dir="{annovar_dir}"
annovar_humandb="{annovar_humandb}"
annovar_protocol="{annovar_protocol}"
annovar_operation="{annovar_operation}"

# tools
fastqc="{fastqc}"
bwa="{bwa}"
samtools="{samtools}"
bcftools="{bcftools}"
java="{java}"
picard="{picard}"
gatk="{gatk}"
vcf_concat="{vcf_concat}"
vcf_sort="{vcf_sort}"
vcf_merge="{vcf_merge}"
bgzip="{bgzip}"
tabix="{tabix}"
