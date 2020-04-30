"""
example:
    #  check run options and generate draft.json
    python wes_pipeline_on_hg19.py

    # it is recommended to rename draft.json to save it from overwriting
    # initial config generation, full version
    python wes_pipeline_on_hg19.py \\
        --draft_settings_json ./__draft_settings__.json \\
        --project_root ./test_alignment \\
        --script_dir_name scripts \\
        --fastq_dirs_list ./test_fastq_gz \\
        --sample_delimiter . \\
        --fastq_extension .fastq.gz \\
        --R1_fastq_extension .R1.fastq.gz \\
        --R2_fastq_extension .R2.fastq.gz \\
        --add_tokens \\
        --debug

    # precise config tuning

    # scripts generation
    python wes_pipeline_on_hg19.py -j /path/to/save/wes/project/scripts/default_settings.json

    # submit scripts to workload manager

    cd /path/to/save/wes/project/scripts/; for i in $( ls *.ss.sh ); do qsub $i; done

    a list of capture targets added to WGS pipeline to extract exomes
    https://gatkforums.broadinstitute.org/gatk/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals
"""

import argparse
import copy
import glob
import json
import os
from collections import defaultdict

__NOT_READY__ = "NOT_READY"
__READY__ = "READY"
__ALMOST_READY__ = "ALMOST_READY"
__DRAFT_SETTINGS_FILE__ = "__draft_settings__.json"


def main():
    settings = parse_arguments_to_settings()
    if settings["ready"] == __READY__:
        run_pipeline(settings)
    elif settings["ready"] == __ALMOST_READY__:
        save_project_settings_json(settings)
    else:
        print(__doc__)
        save_draft_settings_json()
        print("# please provide --project_root and --fastq_dirs_list")


def parse_arguments_to_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", default=None, required=False,
                        help="generated and edited json with all required parameters")
    parser.add_argument("-d", "--draft_settings_json", default=None, required=False,
                        help="path to draft_json with local variables like tools and databases")
    parser.add_argument("-p", "--project_root", default=None, required=False,
                        help="folder to be created that will store pipeline output")
    parser.add_argument("-s", "--script_dir_name", default=None, required=False,
                        help="folder to store all runnable scripts, usually scripts")
    parser.add_argument("-n", "--project_settings_json_name", default=None, required=False,
                        help="folder to be created that will store pipeline output")
    parser.add_argument("-f", "--fastq_dirs_list", default=None, required=False, nargs="+",
                        help="")
    parser.add_argument("--sample_delimiter", default=None, required=False,
                        help="sample delimiter, usually _ (underscore) or . (dot)")
    parser.add_argument("--fastq_extension", default=None, required=False,
                        help="extension of input files, usually .fastq.gz or.fastq")
    parser.add_argument("--R1_fastq_extension", default=None, required=False,
                        help="extension to distinguish R1 from R2")
    parser.add_argument("--R2_fastq_extension", default=None, required=False,
                        help="extension to distinguish R2 from R1")
    parser.add_argument("--add_tokens", action="store_true",
                        help="add tokens to distinguish steps completed, helps on reruns ")
    parser.add_argument("--debug", action="store_true",
                        help="print additional info messages")
    #
    args = parser.parse_args()
    if args.settings_json:
        settings = json.load(open(args.settings_json))  # config exist, import it
        __samples_list__ = settings["samples_list"]  # get list of target samples from config
        __samples_dict__ = load_fastq_samples(settings)  # find all fastq files
        settings["samples_dict"] = {  # filter target samples from all fastq
            list_key: dict_value
            for list_key in __samples_list__
            for dict_key, dict_value in __samples_dict__.items()
            if list_key == dict_key or list_key + "_m" == dict_key
        }
        settings["ready"] = __READY__
    elif args.project_root and args.fastq_dirs_list:
        draft_settings_file = args.draft_settings_json if args.draft_settings_json else __DRAFT_SETTINGS_FILE__
        settings = json.load(open(draft_settings_file))
        project_settings = {
            "draft_settings_json": args.draft_settings_json,
            "project_root": args.project_root,
            "script_dir_name": args.script_dir_name,
            "project_settings_json_name": args.project_settings_json_name,
            "fastq_dirs_list": args.fastq_dirs_list,
            "sample_delimiter": args.sample_delimiter,
            "fastq_extension": args.fastq_extension,
            "R1_fastq_extension": args.R1_fastq_extension,
            "R2_fastq_extension": args.R2_fastq_extension,
            "add_tokens": args.add_tokens,
            "debug": args.debug,
            "ready": __ALMOST_READY__,
        }
        settings.update({k: v for k, v in project_settings.items() if v})
        settings["samples_list"] = sorted(load_fastq_samples(settings))
    else:
        settings = {
            "ready": __NOT_READY__,
        }
    return settings


def load_fastq_samples(settings):
    if settings["debug"]:
        print("# load_fastq_samples: settings\n#", settings)
    #
    fastq_dirs_list = settings["fastq_dirs_list"]
    sample_delimiter = settings["sample_delimiter"]
    fastq_extension = settings["fastq_extension"]
    R1_fastq_extension = settings["R1_fastq_extension"]
    R2_fastq_extension = settings["R2_fastq_extension"]
    #
    sample_dict = defaultdict(lambda: defaultdict(str))
    for fastq_dir in fastq_dirs_list:
        for fastq in glob.iglob(fastq_dir + '/**' + fastq_extension):
            sample = os.path.basename(fastq).split(sample_delimiter)[0]
            if fastq.endswith(R1_fastq_extension):
                sample_dict[sample]["read1"] = fastq
            elif fastq.endswith(R2_fastq_extension):
                sample_dict[sample]["read2"] = fastq
    sample_dict = {
        key: value
        for key, value in sample_dict.items()
        if key + "_m" not in sample_dict
    }  # to exclude unmerged samples
    #
    if settings["debug"]:
        print("# load_fastq_samples: samples\n#", sample_dict)
    return sample_dict


def run_pipeline(settings):
    for sample in sorted(settings["samples_dict"]):
        sample_settings = copy.deepcopy(settings)
        sample_settings["sample"] = sample
        sample_settings = get_settings_for_SSAP(sample_settings)
        cmd_list = get_cmd_list_for_SSAP(sample_settings)
        write_cmd_list_to_file(sample_settings, cmd_list)
    # debug print
    print(f"# ls {settings['project_script_dir']}")
    print(f"# cd {settings['project_script_dir']}")
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do echo $i; done")
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do qsub $i; done")


def save_project_settings_json(settings):
    settings["project_script_dir"] = os.path.join(
        settings["project_root"],
        settings["script_dir_name"],
    )
    mkdir(settings["project_script_dir"])
    project_settings_json_file = os.path.join(
        settings["project_script_dir"],
        settings["project_settings_json_name"],
    )
    with open(project_settings_json_file, "w") as f:
        project_settings_str = json.dumps(settings, indent=4, sort_keys=False)
        f.write(project_settings_str)
    # debug print
    print(f"# ls   {settings['project_script_dir']}")
    print(f"# cd   {settings['project_script_dir']}")
    print(f"# more {project_settings_json_file}")
    print(f"# nano {project_settings_json_file}")
    print(f"# python __script_name__.py -j {project_settings_json_file}")


def save_draft_settings_json():
    def get_draft_settings_json():
        draft_settings_json = {
            "number_of_threads": "4",
            "project_settings_json_name": "project_settings.json",
            "script_dir_name": "scripts",
            "sample_delimiter": "_",
            "fastq_extension": ".fastq.gz",
            "R1_fastq_extension": ".R1.fastq.gz",
            "R2_fastq_extension": ".R2.fastq.gz",

            "add_tokens": False,
            "debug": False,
            "samples_list": [],

            # for picard_ARRG block
            "RGID": "__sample__",
            "RGLB": "__sample__",
            "RGPL": "ILLUMINA",
            "RGPU": "SureSelectV4",
            "RGSM": "__sample__",
            "RGCN": "NLA",

            # tools block
            "fastqc": "fastqc",
            "bwa": "bwa",
            "samtools": "samtools",
            "bcftools": "bcftools",
            "java": "java",
            "picard": "picard",
            "gatk": "gatk",
            "vcf_concat": "vcf-concat",
            "vcf_sort": "vcf-sort",
            "vcf_merge:": "vcf-merge",
            "bgzip": "bgzip",
            "tabix": "tabix",

            # databases block
            "ref": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
            "ref_ucsc_hg19": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
            "ref_mtDNA": "/home/PublicData/h.sapiens_mtDNA/HS_mtDNA.fa",
            "gold_indel": "/home/PublicData/broadinstitute/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
            "oneKG_indel": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf",
            "oneKG_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf",
            "dbsnp": "/home/PublicData/broadinstitute/2.8/hg19/dbsnp_138.hg19.vcf",
            "onmi_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_omni2.5.hg19.sites.vcf",
            "hapmap_snp": "/home/PublicData/broadinstitute/2.8/hg19/hapmap_3.3.hg19.sites.vcf",
            "target_region": "/home/PublicData/Agilent_v4_71m_reduced.bed",
        }
        return draft_settings_json

    draft_settings_json = get_draft_settings_json()
    with open(__DRAFT_SETTINGS_FILE__, "w") as f:
        draft_json_str = json.dumps(draft_settings_json, indent=4, sort_keys=False)
        f.write(draft_json_str)
    print(f"# more {__DRAFT_SETTINGS_FILE__}")
    print(f"# nano {__DRAFT_SETTINGS_FILE__}")


def mkdir(dir_name):
    os.makedirs(dir_name, exist_ok=True)


def get_settings_for_SSAP(sample_dict):
    sample = sample_dict["sample"]
    sample_dir = os.path.join(sample_dict["project_root"], sample)
    sample_path_prefix = os.path.join(sample_dir, sample)

    _dict = {
        "sample": sample,
        "sample_dir": sample_dir,

        "read1": sample_dict["samples_dict"][sample]["read1"],
        "read2": sample_dict["samples_dict"][sample]["read2"],

        "sam_bwa_mem": sample_path_prefix + ".bwa_mem.sam",
        "bam_samtools_view": sample_path_prefix + ".samtools_view.bam",
        "bam_samtools_sort": sample_path_prefix + ".samtools_sort.bam",
        "bam_picard_ARRG": sample_path_prefix + ".picard_ARRG.bam",

        "RGID": sample if sample_dict["RGID"] == "__sample__" else sample_dict["RGID"],
        "RGLB": sample if sample_dict["RGLB"] == "__sample__" else sample_dict["RGLB"],
        "RGPL": sample if sample_dict["RGPL"] == "__sample__" else sample_dict["RGPL"],
        "RGPU": sample if sample_dict["RGPU"] == "__sample__" else sample_dict["RGPU"],
        "RGSM": sample if sample_dict["RGSM"] == "__sample__" else sample_dict["RGSM"],
        "RGCN": sample if sample_dict["RGCN"] == "__sample__" else sample_dict["RGCN"],

        "bam_picard_MD": sample_path_prefix + ".picard_MD.bam",
        "txt_picard_MD": sample_path_prefix + ".picard_MD.metrix.txt",

        "bam_picard_SS_SNAUT": sample_path_prefix + ".picard_SS_SNAUT.bam",

        "table_gatk_BR": sample_path_prefix + ".gatk_BR.table",
        "table_gatk_BR_BQSR": sample_path_prefix + ".gatk_BR_BQSR.table",
        "pdf_gatk_AC": sample_path_prefix + ".gatk_AC.pdf",
        "bam_gatk_PR_BR_BQSR": sample_path_prefix + ".gatk_PR_BR_BQSR.bam",

        "vcf_gatk_HC": sample_path_prefix + ".gatk_HC.vcf",
        "bam_gatk_HC": sample_path_prefix + ".gatk_HC.bam",

        "recal_gatk_VR_SNP": sample_path_prefix + ".gatk_VR_SNP.recal",
        "tranches_gatk_VR_SNP": sample_path_prefix + ".gatk_VR_SNP.tranches",
        "plots_gatk_VR_SNP": sample_path_prefix + ".gatk_VR_SNP_plots.Rscript",
        "vcf_gatk_AR_SNP": sample_path_prefix + ".gatk_AR_SNP.vcf",

        "recal_gatk_VR_INDEL": sample_path_prefix + ".gatk_VR_INDEL.recal",
        "tranches_gatk_VR_INDEL": sample_path_prefix + ".gatk_VR_INDEL.tranches",
        "plots_gatk_VR_INDEL": sample_path_prefix + ".gatk_VR_INDEL_plots.Rscript",
        "vcf_gatk_AR_INDEL": sample_path_prefix + ".gatk_AR_INDEL.vcf",

        "vcf_gatk_SV_SNP_raw": sample_path_prefix + ".gatk_SV_SNP_raw.vcf",
        "vcf_gatk_SV_SNP_fil": sample_path_prefix + ".gatk_SV_SNP_fil.vcf",

        "vcf_gatk_SV_INDEL_raw": sample_path_prefix + ".gatk_SV_INDEL_raw.vcf",
        "vcf_gatk_SV_INDEL_fil": sample_path_prefix + ".gatk_SV_INDEL_fil.vcf",

        "vcf_vcftools_concat": sample_path_prefix + ".vcftools_concat.FINAL.vcf",
        "vcf_vcftools_sorted": sample_path_prefix + ".vcftools_sorted.FINAL.vcf",
    }
    sample_dict.update(_dict)
    return sample_dict


###############################################################################
def get_cmd_list_for_SSAP(sample_settings):
    if sample_settings["debug"]:
        print("\n\nget_cmd_list_for_SSAP", sample_settings)
    cmd_list = [
        get_mkdir_cmd(sample_settings),
        get_cmd_bwa_mem_sam(sample_settings),
        get_cmd_samtools_view_bam(sample_settings),
        get_cmd_samtools_sort_bam(sample_settings),
        get_cmd_picard_ARRG_bam(sample_settings),
        get_cmd_picard_MD_bam(sample_settings),

        get_cmd_gatk_BR_table(sample_settings),
        get_cmd_gatk_BR_BQSR_table(sample_settings),
        get_cmd_gatk_AC_pdf(sample_settings),
        get_cmd_gatk_PR_bam(sample_settings),

        get_cmd_gatk_HC_vcf(sample_settings),

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
def get_mkdir_cmd(d):
    return "mkdir -p {sample_dir}".format(**d)


###############################################################################
def reduce_spaces_and_newlines(func):
    def wrapper(*args, **kwargs):
        s = func(*args, **kwargs)
        s = s.replace("\n", " ")
        s = " ".join([i for i in s.split(" ") if i])
        return s

    return wrapper


@reduce_spaces_and_newlines
def get_cmd(d):
    if d["add_tokens"]:
        d["token"] = "{sample_dir}/token.{sample}.{token_suffix}".format(**d)
        d["flags"] = " && ".join([" [ -f {:s} ] ".format(i) for i in d["files_list"]]) + \
                     " && [ ! -f {token} ] ".format(**d)
        cmd = """
            {flags} &&
            dt1=`date +%y%m%d_%H%M%S` && 
            echo $dt1 {token} &&
            {main_cmd} &&
            du {out_file} > {out_file}.$dt1.du &&
            md5sum {out_file} > {out_file}.$dt1.md5 &&
            dt2=`date +%y%m%d_%H%M%S` &&
            echo $dt1 $dt2 > {token} ||
            echo "TOKEN SKIPPED {token}"
            """.format(**d)
    else:
        cmd = "{main_cmd}".format(**d)
    return cmd


@reduce_spaces_and_newlines
def clear_after_competion(d):
    d["token_suffix"] = "vcftools_sort_vcf"
    d["token"] = "{sample_dir}/token.{sample}.{token_suffix}".format(**d)
    if d["add_tokens"]:
        cmd = """ [ -f {token} ] && [ -f {vcftools_sorted_vcf} ] &&
            rm -f
            {sam} {bam} {sorted_bam} {ARRG_bam}
            {IR_bam}
            {gatk_AR_SNP_vcf} {gatk_AR_INDEL_vcf}
            {gatk_SV_SNP_raw_vcf}  {gatk_SV_INDEL_raw_vcf}
            {vcftools_concat_vcf}
            """.format(**d)
        # {gatk_SV_SNP_fil_vcf} {gatk_SV_INDEL_fil_vcf} # for vcftools concatenate - sometimes library doesnt loaded well
    else:
        cmd = ""
    return cmd


###############################################################################
def get_cmd_bwa_mem_sam(d):
    d["token_suffix"] = "read1_read2_2_sam_bwa_mem"
    d["files_list"] = [d["read1"], d["read2"]]
    d["out_file"] = d["sam_bwa_mem"]
    d["main_cmd"] = bash_bwa_mem_sam(d)
    return get_cmd(d)


def bash_bwa_mem_sam(d):
    return """{bwa} mem -M -t {number_of_threads} {ref} {read1} {read2} > {sam_bwa_mem}""".format(**d)


###############################################################################
def get_cmd_samtools_view_bam(d):
    d["token_suffix"] = "sam_2_bam_samtools_view"
    d["files_list"] = [d["sam_bwa_mem"]]
    d["out_file"] = d["bam_samtools_view"]
    d["main_cmd"] = bash_samtools_view_bam(d)
    return get_cmd(d)


def bash_samtools_view_bam(d):
    return """{samtools} view -bT {ref} {sam_bwa_mem} > {bam_samtools_view}""".format(**d)


###############################################################################
def get_cmd_samtools_sort_bam(d):
    d["token_suffix"] = "bam_2_bam_samtools_sort"
    d["files_list"] = [d["bam_samtools_view"]]
    d["out_file"] = d["bam_samtools_sort"]
    d["main_cmd"] = bash_samtools_sort_bam(d)
    return get_cmd(d)


def bash_samtools_sort_bam(d):
    return """{samtools} sort -l 9 -O bam {bam_samtools_view} > {bam_samtools_sort}""".format(**d)


###############################################################################
def get_cmd_picard_ARRG_bam(d):
    d["token_suffix"] = "bam_2_bam_picard_ARRG"
    d["files_list"] = [d["bam_samtools_sort"]]
    d["out_file"] = d["bam_picard_ARRG"]
    d["main_cmd"] = bash_picard_ARRG_bam(d)
    return get_cmd(d)


def bash_picard_ARRG_bam(d):
    return """
        {picard}
        AddOrReplaceReadGroups
        INPUT={bam_samtools_sort}
        OUTPUT={bam_picard_ARRG}
        SORT_ORDER=coordinate
        RGID={RGID}
        RGLB={RGLB}
        RGPL={RGPL}
        RGPU={RGPU}
        RGSM={RGSM}
        RGCN={RGCN}
        CREATE_INDEX=true
        VALIDATION_STRINGENCY=LENIENT
        MAX_RECORDS_IN_RAM=1000000
        """.format(**d)
    # TMP_DIR={tmp_dir}


###############################################################################
def get_cmd_picard_MD_bam(d):
    d["token_suffix"] = "picard_ARRG_bam_2_MD_bam"
    d["files_list"] = [d["bam_picard_ARRG"]]
    d["out_file"] = d["bam_picard_MD"]
    d["main_cmd"] = bash_picard_MD_bam(d)
    return get_cmd(d)


def bash_picard_MD_bam(d):
    return """
        {picard}
        MarkDuplicates
        INPUT={bam_picard_ARRG}
        OUTPUT={bam_picard_MD}
        METRICS_FILE={txt_picard_MD}
        ASSUME_SORTED=true
        CREATE_INDEX=true
        VALIDATION_STRINGENCY=LENIENT
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
        """.format(**d)
    # TMP_DIR={tmp_dir}


###############################################################################
def get_cmd_gatk_RTC_intervals(d):
    pass
    # https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847


def bash_gatk_RTC_intervals(d):
    pass
    # https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847


###############################################################################
def get_cmd_gatk_IR_bam(d):
    pass
    # https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847


def bash_gatk_IR_bam(d):
    pass
    # https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847


###############################################################################
def get_cmd_gatk_BR_table(d):
    d["token_suffix"] = "bam_2_table_gatk_BR"
    d["files_list"] = [d["bam_picard_MD"]]
    d["out_file"] = d["table_gatk_BR"]
    d["main_cmd"] = bash_gatk_BR_table(d)
    return get_cmd(d)


def bash_gatk_BR_table(d):
    # https://gatkforums.broadinstitute.org/gatk/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals
    return """
        {gatk}
        -T BaseRecalibrator
        -R {ref}
        -I {bam_picard_MD}
        -L {target_region}
        -knownSites {dbsnp}
        -knownSites {gold_indel}
        -knownSites {oneKG_indel}
        -o {table_gatk_BR}
        """.format(**d)


###############################################################################
def get_cmd_gatk_BR_BQSR_table(d):
    d["token_suffix"] = "bam_table_2_table_gatk_BR_BQSR"
    d["files_list"] = [d["bam_picard_MD"], d["table_gatk_BR"]]
    d["out_file"] = d["table_gatk_BR_BQSR"]
    d["main_cmd"] = bash_gatk_BR_BQSR_table(d)
    return get_cmd(d)


def bash_gatk_BR_BQSR_table(d):
    # https://gatkforums.broadinstitute.org/gatk/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals
    return """
        {gatk}
        -T BaseRecalibrator
        -R {ref}
        -I {bam_picard_MD}
        -L {target_region}
        -knownSites {dbsnp}
        -knownSites {gold_indel}
        -knownSites {oneKG_indel}
        -BQSR {table_gatk_BR}
        -o {table_gatk_BR_BQSR}
        """.format(**d)


###############################################################################
def get_cmd_gatk_AC_pdf(d):
    d["token_suffix"] = "table_table_2_pdf_gatk_AC"
    d["files_list"] = [d["table_gatk_BR"], d["table_gatk_BR_BQSR"]]
    d["out_file"] = d["pdf_gatk_AC"]
    d["main_cmd"] = bash_gatk_AC_pdf(d)
    return get_cmd(d)


def bash_gatk_AC_pdf(d):
    return """
        {gatk}
        -T AnalyzeCovariates
        -R {ref}
        -before {table_gatk_BR}
        -after {table_gatk_BR_BQSR}
        -plots {pdf_gatk_AC}
        """.format(**d)


###############################################################################
def get_cmd_gatk_PR_bam(d):
    d["token_suffix"] = "bam_table_2_bam_gatk_PR"
    d["files_list"] = [d["bam_picard_MD"], d["table_gatk_BR_BQSR"]]
    d["out_file"] = d["bam_gatk_PR_BR_BQSR"]
    d["main_cmd"] = bash_gatk_PR_bam(d)
    return get_cmd(d)


def bash_gatk_PR_bam(d):
    return """
        {gatk}
        -T PrintReads
        -R {ref}
        -I {bam_picard_MD}
        -BQSR {table_gatk_BR_BQSR}
        -o {bam_gatk_PR_BR_BQSR}
        """.format(**d)


###############################################################################
def get_cmd_gatk_HC_vcf(d):
    d["token_suffix"] = "bam_2_vcf_gatk_HC"
    d["files_list"] = [d["bam_gatk_PR_BR_BQSR"]]
    d["out_file"] = d["vcf_gatk_HC"]
    d["main_cmd"] = bash_gatk_HC_vcf(d)
    return get_cmd(d)


def bash_gatk_HC_vcf(d):
    # https://gatkforums.broadinstitute.org/gatk/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals
    return """
        {gatk}
        -T HaplotypeCaller
        -R {ref}
        -I {bam_gatk_PR_BR_BQSR}
        -D {dbsnp}
        -L {target_region}
        --genotyping_mode DISCOVERY
        -stand_call_conf 30
        -o {vcf_gatk_HC}
        -nct {number_of_threads}
        """.format(**d)


###############################################################################
def get_cmd_gatk_UG_vcf(d):  # not tested, not ready
    d["token_suffix"] = "bam_2_vcf_gatk_UG"
    d["files_list"] = [d["bam_gatk_PR_BR_BQSR"]]
    d["out_file"] = d["vcf_gatk_UG"]
    d["main_cmd"] = bash_gatk_UG_vcf(d)
    return get_cmd(d)


def bash_gatk_UG_vcf(d):  # not tested, not ready
    return """
        {gatk}
        -T UnifiedGenotyper
        -R {ref}
        -I {bam_gatk_PR_BR_BQSR}
        -D {dbsnp}
        --genotyping_mode DISCOVERY
        -o {vcf_gatk_UG}
        --output_mode EMIT_ALL_SITES
        -nct {number_of_threads}
        --max_num_PL_values 100
        --min_base_quality_score 17
        """.format(**d)


###############################################################################
def get_cmd_gatk_VR_SNP_recal_tranches(d):
    d["token_suffix"] = "vcf_2_recal_tranches_gatk_VR"
    d["files_list"] = [d["vcf_gatk_HC"]]
    d["out_file"] = d["recal_gatk_VR_SNP"]
    d["main_cmd"] = bash_gatk_VR_SNP_recal_tranches(d)
    return get_cmd(d)


def bash_gatk_VR_SNP_recal_tranches(d):
    return """
        {gatk}
        -T VariantRecalibrator
        -R {ref}
        -input {vcf_gatk_HC}
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap_snp}
        -resource:omni,known=false,training=true,truth=true,prior=12.0 {onmi_snp}
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 {oneKG_snp}
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp}
        -an DP
        -an QD
        -an FS
        -an MQRankSum
        -an ReadPosRankSum
        -mode SNP
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
        --maxGaussians 4
        -recalFile {recal_gatk_VR_SNP}
        -tranchesFile {tranches_gatk_VR_SNP}
        -rscriptFile {plots_gatk_VR_SNP}
        """.format(**d)


###############################################################################
def get_cmd_gatk_AR_SNP_vcf(d):
    d["token_suffix"] = "vcf_recal_tranches_2_vcf_gatk_AR_SNP"
    d["files_list"] = [d["vcf_gatk_HC"], d["recal_gatk_VR_SNP"], d["tranches_gatk_VR_SNP"]]
    d["out_file"] = d["vcf_gatk_AR_SNP"]
    d["main_cmd"] = bash_gatk_AR_SNP_vcf(d)
    return get_cmd(d)


def bash_gatk_AR_SNP_vcf(d):
    return """
        {gatk}
        -T ApplyRecalibration
        -R {ref}
        -input {vcf_gatk_HC}
        -mode SNP
        --ts_filter_level 99.0
        -recalFile {recal_gatk_VR_SNP}
        -tranchesFile {tranches_gatk_VR_SNP}
        -o {vcf_gatk_AR_SNP}
        """.format(**d)


###############################################################################
def get_cmd_gatk_VR_INDEL_recal_tranches(d):
    d["token_suffix"] = "vcf_2_recal_tranches_gatk_VR_INDEL"
    d["files_list"] = [d["vcf_gatk_AR_SNP"]]
    d["out_file"] = d["recal_gatk_VR_INDEL"]
    d["main_cmd"] = bash_gatk_VR_INDEL_recal_tranches(d)
    return get_cmd(d)


def bash_gatk_VR_INDEL_recal_tranches(d):
    return """
        {gatk}
        -T VariantRecalibrator
        -R {ref}
        -input {vcf_gatk_AR_SNP}
        -resource:mills,known=true,training=true,truth=true,prior=12.0 {gold_indel}
        -an DP
        -an FS
        -an MQRankSum
        -an ReadPosRankSum
        -mode INDEL
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
        --maxGaussians 4
        -recalFile {recal_gatk_VR_INDEL}
        -tranchesFile {tranches_gatk_VR_INDEL}
        -rscriptFile {plots_gatk_VR_INDEL}
        """.format(**d)


###############################################################################
def get_cmd_gatk_AR_INDEL_vcf(d):
    d["token_suffix"] = "vcf_recal_tranches_2_vcf_gatk_AR_INDEL"
    d["files_list"] = [d["vcf_gatk_AR_SNP"], d["recal_gatk_VR_INDEL"], d["tranches_gatk_VR_INDEL"]]
    d["out_file"] = d["vcf_gatk_AR_INDEL"]
    d["main_cmd"] = bash_gatk_AR_INDEL_vcf(d)
    return get_cmd(d)


def bash_gatk_AR_INDEL_vcf(d):
    return """
        {gatk}
        -T ApplyRecalibration
        -R {ref}
        -input {vcf_gatk_AR_SNP}
        -mode INDEL
        --ts_filter_level 99.0
        -recalFile {recal_gatk_VR_INDEL}
        -tranchesFile {tranches_gatk_VR_INDEL}
        -o {vcf_gatk_AR_INDEL}
        """.format(**d)


###############################################################################
def get_cmd_gatk_SV_SNP_raw_vcf(d):
    d["token_suffix"] = "vcf_2_vcf_gatk_SV_raw"
    d["files_list"] = [d["vcf_gatk_AR_INDEL"]]
    d["out_file"] = d["vcf_gatk_SV_SNP_raw"]
    d["main_cmd"] = bash_gatk_SV_SNP_raw_vcf(d)
    return get_cmd(d)


def bash_gatk_SV_SNP_raw_vcf(d):
    return """
        {gatk}
        -T SelectVariants
        -R {ref}
        -V {vcf_gatk_AR_INDEL}
        -selectType SNP
        -o {vcf_gatk_SV_SNP_raw}
        """.format(**d)


###############################################################################
def get_cmd_gatk_VF_SNP_fil_vcf(d):
    d["token_suffix"] = "vcf_2_vcf_gatk_SV_SNP"
    d["files_list"] = [d["vcf_gatk_SV_SNP_raw"]]
    d["out_file"] = d["vcf_gatk_SV_SNP_fil"]
    d["main_cmd"] = bash_gatk_VF_SNP_fil_vcf(d)
    return get_cmd(d)


def bash_gatk_VF_SNP_fil_vcf(d):
    return """
        {gatk}
        -T VariantFiltration
        -R {ref}
        -V {vcf_gatk_SV_SNP_fil}
        --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        --filterName "SNP_FAIL"
        -o {vcf_gatk_SV_SNP_fil}
        """.format(**d)


###############################################################################
def get_cmd_gatk_SV_INDEL_raw_vcf(d):
    d["token_suffix"] = "vcf_2_vcf_SV_INDEL_raw"
    d["files_list"] = [d["vcf_gatk_AR_INDEL"]]
    d["out_file"] = d["vcf_gatk_SV_INDEL_raw"]
    d["main_cmd"] = bash_gatk_SV_INDEL_raw_vcf(d)
    return get_cmd(d)


def bash_gatk_SV_INDEL_raw_vcf(d):
    return """
        {gatk}
        -T SelectVariants
        -R {ref}
        -V {vcf_gatk_AR_INDEL}
        -selectType INDEL
        -o {vcf_gatk_SV_INDEL_raw}
        """.format(**d)


###############################################################################
def get_cmd_gatk_VF_INDEL_fil_vcf(d):
    d["token_suffix"] = "vcf_2_vcf_gatk_VF_INDEL"
    d["files_list"] = [d["vcf_gatk_SV_INDEL_raw"]]
    d["out_file"] = d["vcf_gatk_SV_INDEL_fil"]
    d["main_cmd"] = bash_gatk_VF_INDEL_fil_vcf(d)
    return get_cmd(d)


def bash_gatk_VF_INDEL_fil_vcf(d):
    return """
        {gatk}
        -T VariantFiltration
        -R {ref}
        -V {vcf_gatk_SV_INDEL_raw}
        --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        --filterName "INDEL_FAIL"
        -o {vcf_gatk_SV_INDEL_fil}
        """.format(**d)


###############################################################################
def get_cmd_vcf_concat(d):
    d["token_suffix"] = "vcf_vcf_2_vcf_vcftools_concat"
    d["files_list"] = [d["vcf_gatk_SV_SNP_fil"], d["vcf_gatk_SV_INDEL_fil"]]
    d["out_file"] = d["vcf_vcftools_concat"]
    d["main_cmd"] = bash_vcf_concat(d)
    return get_cmd(d)


def bash_vcf_concat(d):
    return """{vcf_concat} {vcf_gatk_SV_SNP_fil} {vcf_gatk_SV_INDEL_fil} > {vcf_vcftools_concat}""".format(**d)


###############################################################################
def get_cmd_vcf_sort(d):
    d["token_suffix"] = "vcf_2_vcf_vcftools_sort"
    d["files_list"] = [d["vcf_vcftools_concat"]]
    d["out_file"] = d["vcf_vcftools_sorted"]
    d["main_cmd"] = bash_vcf_sort(d)
    return get_cmd(d)


def bash_vcf_sort(d):
    return """cat {vcf_vcftools_concat} | {vcf_sort} > {vcf_vcftools_sorted}""".format(**d)


###############################################################################
def write_cmd_list_to_file(sample_settings, cmd_list):
    script_file = os.path.join(
        sample_settings["project_script_dir"],
        sample_settings["sample"] + ".ss.sh",
    )
    print("###")
    with open(script_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for _, cmd in enumerate(cmd_list):
            print("#", _, cmd)
            new_line = cmd + "\n\n"
            f.write(new_line)


###############################################################################
if __name__ == "__main__":
    main()
