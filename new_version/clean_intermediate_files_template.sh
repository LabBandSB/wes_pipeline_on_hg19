###################################################################
# delete bwa_mem
###################################################################
token1="${alignment_dir}/token.${sample}.fastq_2_sam_bwa_mem"
token2="${alignment_dir}/token.${sample}.sam_2_bam_samtools_view"
output_file="${alignment_dir}/${sample}.bwa_mem.sam"
[ -f ${token1} ] && \
[ -f ${token2} ] && \
[ -f ${output_file} ] && \
rm -f ${output_file}

###################################################################
# delete samtools_view
###################################################################
token1="${alignment_dir}/token.${sample}.sam_2_bam_samtools_view"
token2="${alignment_dir}/token.${sample}.bam_2_bam_samtools_sort"
output_file="${alignment_dir}/${sample}.samtools_view.bam"
[ -f ${token1} ] && \
[ -f ${token2} ] && \
[ -f ${output_file} ] && \
rm -f ${output_file}

###################################################################
# delete samtools_sort
###################################################################
# token1="${alignment_dir}/token.${sample}.bam_2_bam_samtools_sort"
# token2="${alignment_dir}/token.${sample}.bam_2_bam_picard_ARRG"
# output_file="${alignment_dir}/${sample}.samtools_sort.bam"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

###################################################################
# delete picard_ARRG
###################################################################
token1="${alignment_dir}/token.${sample}.bam_2_bam_picard_ARRG"
token2="${alignment_dir}/token.${sample}.bam_2_bam_picard_MD"
output_file="${alignment_dir}/${sample}.picard_ARRG.bam"
[ -f ${token1} ] && \
[ -f ${token2} ] && \
[ -f ${output_file} ] && \
rm -f ${output_file}

# ###################################################################
# # delete BR_table
# ###################################################################
# token1="${alignment_dir}/token.${sample}.bam_2_bam_picard_MD"
# token2="${alignment_dir}/token.${sample}.bam_2_txt_gatk_BR"
# output_file="${alignment_dir}/${sample}.BR_table.txt"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# ###################################################################
# # delete picard_ARRG
# ###################################################################
# token1="${alignment_dir}/token.${sample}.bam_2_txt_gatk_BR"
# token2="${alignment_dir}/token.${sample}.txt_2_txt_gatk_BR"
# output_file="${alignment_dir}/${sample}.BQSR_BR_table.txt"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

###################################################################
# delete gatk_PR_BQSR_BR_table
###################################################################
token1="${alignment_dir}/token.${sample}.bam_2_bam_gatk_PR"
token2="${alignment_dir}/token.${sample}.bam_2_vcf_gatk_HC"
output_file="${alignment_dir}/${sample}.gatk_PR_BQSR_BR_table.bam"
[ -f ${token1} ] && \
[ -f ${token2} ] && \
[ -f ${output_file} ] && \
rm -f ${output_file}

# ###################################################################
# # delete gatk_HC
# ###################################################################
# token1="${alignment_dir}/token.${sample}.bam_2_vcf_gatk_HC"
# token2="${alignment_dir}/token.${sample}.vcf_2_FINAL"
# output_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# output_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.gatk_AR_INDEL_VQSR.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# output_file="${alignment_dir}/${sample}.gatk_SV_SNP.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# output_file="${alignment_dir}/${sample}.gatk_VF_SNP.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# output_file="${alignment_dir}/${sample}.gatk_SV_INDEL.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# output_file="${alignment_dir}/${sample}.gatk_VF_INDEL.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# output_file="${alignment_dir}/${sample}.vcf_concat.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}

# output_file="${alignment_dir}/${sample}.vcf_sort.vcf"
# [ -f ${token1} ] && \
# [ -f ${token2} ] && \
# [ -f ${output_file} ] && \
# rm -f ${output_file}
