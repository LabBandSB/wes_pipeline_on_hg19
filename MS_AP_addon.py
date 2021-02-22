# "project_root": "/home/adminrig/projects/20210208_onco_patient/alignment"
def run_MS_pipeline(settings):
    input_bams = ''
    for bam in get_files_generator([settings['project_root'], ".BQSR_BR.bam"]):
        input_bams += ' -I '+ bam

    sample_settings = copy.deepcopy(settings)
    sample_settings["sample"] = sample_settings['gatk_MSHC_name']
    sample_settings = get_settings_for_SSAP(sample_settings)
    sample_settings['input_bams'] = input_bams
    cmd_list = get_cmd_list_for_MSAP(sample_settings)
    write_cmd_list_to_file(sample_settings, cmd_list)

def get_files_generator(dirs_list, extension):
        for _dir in sorted(dirs_list):
            for _file in sorted(os.listdir(_dir)):
                if not _file:
                    continue
                _file_path = os.path.join(_dir, _file)
                if os.path.isfile(_file_path) and _file_path.endswith(extension):
                    yield _file_path                
                elif os.path.isdir(_file_path):
                    yield from get_files_generator([_file_path], extension)
    

def get_cmd_list_for_MSAP(sample_settings):
    cmd_list = [
        get_mkdir_cmd(sample_settings),

        get_cmd_gatk_MSHC_vcf(sample_settings),

        get_cmd_gatk_VR_SNP_recal_tranches(sample_settings),
        get_cmd_gatk_AR_SNP_vcf(sample_settings),

        get_cmd_gatk_VR_INDEL_recal_tranches(sample_settings),
        get_cmd_gatk_AR_INDEL_vcf(sample_settings),

        get_cmd_gatk_SV_SNP_raw_vcf(sample_settings),
        get_cmd_gatk_VF_SNP_fil_vcf(sample_settings),

        get_cmd_gatk_SV_INDEL_raw_vcf(sample_settings),
        get_cmd_gatk_VF_INDEL_fil_vcf(sample_settings),

        get_cmd_vcf_concat(sample_settings),
        get_cmd_vcf_sort(sample_settings),
    ]
    return cmd_list

###############################################################################
def get_cmd_gatk_MSHC_vcf(d):
    d["token_suffix"] = "bam_2_vcf_gatk_HC"
    d["files_list"] = []
    d["out_file"] = d["vcf_gatk_HC"]
    d["main_cmd"] = bash_gatk_MSHC_vcf(d)
    return get_cmd(d)


def bash_gatk_MSHC_vcf(d):
    # https://gatkforums.broadinstitute.org/gatk/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals
    return """
        {gatk}
        -T HaplotypeCaller
        -R {ref}
        {input_bams}
        --dbsnp {dbsnp}
        -L {target_region}
        --genotyping_mode DISCOVERY
        -stand_call_conf 30
        -o {vcf_gatk_HC}
        -nct {number_of_threads}
        """.format(**d)
