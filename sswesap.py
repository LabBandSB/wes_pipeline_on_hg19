"""
example:
    #  check run options
    python sswesap.py


    # initial config generation
    python sswesap.py \\
        --project_root ~/projects/KAZ_WE/KAZ_WE_hg19 \\
        --fastq_dirs_list ~/icebox/fastq_gz/KAZ_WE/ \\
        --sample_delimiter . \\
        --fastq_extension .fastq.gz \\
        --R1_fastq_extension .R1.fastq.gz \\
        --R2_fastq_extension .R2.fastq.gz \\
        --script_dir_name scripts

    # optional: generate config with --add_tokens option
    # to add tokens to each step, for rerun from last failed step


    # precise config tuning

    # scripts generation
    python sswesap.py -j ~/projects/KAZ_WE/KAZ_WE_hg19/scripts/default_settings.json

    # submit scripts to workload manager

    cd ~/projects/KAZ_WE/KAZ_WE_hg19/scripts/; for i in $( ls *.ss.sh ); do qsub $i; done

    # optional: generate config with --run_annovar key

"""
__VERSION__ = "0.0.4"

import os
import argparse
import json

from collections import defaultdict


__NOT_READY__ = "NOT_READY"
__READY__ = "READY"
__ALMOST_READY__ = "ALMOST_READY"


def main():
    settings = parse_arguments_to_settings()
    if settings["ready"] == __ALMOST_READY__:
        save_settings(settings)
    elif settings["ready"] == __READY__:
        run_pipeline(settings)
    else:
        print(__doc__)


def parse_arguments_to_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", default=None, required=False)
    parser.add_argument("--project_root", default=None, required=False)
    parser.add_argument("--fastq_dirs_list", default=[], required=False, nargs="+")
    parser.add_argument("--sample_delimiter", default="_", required=False)
    parser.add_argument("--fastq_extension", default=".fastq.gz", required=False)
    parser.add_argument("--R1_fastq_extension", default=".R1.fastq.gz", required=False)
    parser.add_argument("--R2_fastq_extension", default=".R2.fastq.gz", required=False)
    parser.add_argument("--script_dir_name", default="scripts", required=False)
    parser.add_argument("--run_annovar", action="store_true")
    parser.add_argument("--add_tokens", action="store_true")
    parser.add_argument("--debug", action="store_true")
    #
    args = parser.parse_args()
    if args.settings_json:
        settings = json.load(open(args.settings_json))   # config exist, import it
        __samples_dict__ = load_fastq_samples(settings)  # find all fastq files
        __samples_list__ = settings["samples_list"]  # get list of target samples from config
        settings["samples_dict"] = {  # filter target samples from all fastq
            list_key: dict_value
            for list_key in __samples_list__
            for dict_key, dict_value in __samples_dict__.items()
            if list_key == dict_key or list_key + "_m" == dict_key
        }
        settings["ready"] = __READY__
    elif args.project_root:
        settings = {
            "settings_json": args.settings_json,
            "project_root": args.project_root,
            "fastq_dirs_list": args.fastq_dirs_list,
            "sample_delimiter": args.sample_delimiter,
            "fastq_extension": args.fastq_extension,
            "R1_fastq_extension": args.R1_fastq_extension,
            "R2_fastq_extension": args.R2_fastq_extension,
            "script_dir_name": args.script_dir_name,
            "run_annovar": args.run_annovar,
            "add_tokens": args.add_tokens,
            "debug": args.debug,
            "ready":__ALMOST_READY__,
        }
    else:
        settings = {
            "ready":__NOT_READY__,
        }
    return settings


def load_fastq_samples(settings):
    fastq_dirs_list = settings["fastq_dirs_list"]
    sample_delimiter = settings["sample_delimiter"]
    fastq_extension = settings["fastq_extension"]
    R1_fastq_extension = settings["R1_fastq_extension"]
    R2_fastq_extension = settings["R2_fastq_extension"]
    #
    res = defaultdict(lambda: defaultdict(str))
    for fastq in get_files_generator(fastq_dirs_list, fastq_extension):
        sample = os.path.basename(fastq).split(sample_delimiter)[0]
        if fastq.endswith(R1_fastq_extension):
            res[sample]["read1"] = fastq
        elif fastq.endswith(R2_fastq_extension):
            res[sample]["read2"] = fastq
    res = {
        key: value
        for key, value in res.items()
        if key + "_m" not in res
    }
    return res


def get_files_generator(dirs_list, extension=""):
    for path in dirs_list:
        for data_file in os.listdir(path):
            if data_file:
                data_path = os.path.join(path, data_file)
                if os.path.isfile(data_path) and data_path.endswith(extension):
                    yield data_path
                elif os.path.isdir(data_path):
                    yield from get_files_generator([data_path], extension)


def save_settings(settings):
    settings["project_script_dir"] = os.path.join(
        settings["project_root"],
        settings["script_dir_name"],
    )
    settings.update(get_default_settings(settings))
    mkdir(settings["project_script_dir"])
    json_file = os.path.join(
        settings["project_script_dir"],
        "default_settings.json",
    )
    with open(json_file, "w") as f:
        default_settings_str = json.dumps(settings, indent=4, sort_keys=True)
        f.write(default_settings_str)
    # debug print
    print(f"# ls  {settings['project_script_dir']}")
    print(f"# more {json_file}")
    print(f"# nano {json_file}")
    print(f"# python sswesap.py -j {json_file}")


def get_default_settings(d):
    default_settings_dict = {
        "number_of_threads": "4",
        "fastq_dirs_list": d["fastq_dirs_list"],
        "sample_delimiter": d["sample_delimiter"],
        "fastq_extension": d["fastq_extension"],
        "R1_fastq_extension": d["R1_fastq_extension"],
        "R2_fastq_extension": d["R2_fastq_extension"],
        "samples_list": sorted(load_fastq_samples(d)) if  d["fastq_dirs_list"] else [],
        "project_root": d["project_root"],
        "read_group": {
            "RGID": "__sample__",
            "RGLB": "__sample__",
            "RGPL": "ILLUMINA",
            "RGPU": "SureSelectV4",
            "RGSM": "__sample__",
            "RGCN": "NLA"
        },
        "tools": {
            "fastqc": "",
            "bwa": "/home/Pipeline/bwa/bwa-0.7.12/bwa",
            "samtools": "/home/Pipeline/samtools/samtools-1.2/samtools",
            "bcftools": "/home/Pipeline/bcftools/bcftools-1.2/bcftools",
            "java": "/home/Pipeline/java/jdk1.8.0_101/bin/java",
            "picard": "/home/Pipeline/java/jdk1.8.0_101/bin/java -jar /home/Pipeline/picard/picard_1_130/picard.jar",
            "gatk": "/home/Pipeline/java/jdk1.8.0_101/bin/java -jar /home/Pipeline/gatk/GATK_3_8_1/GenomeAnalysisTK.jar",
            "vcf_concat": "/home/Pipeline/vcftools/vcftools_0.1.13/bin/vcf-concat",
            "vcf_sort": "/home/Pipeline/vcftools/vcftools_0.1.13/bin/vcf-sort",
            "vcf_merge:": "/home/Pipeline/vcftools/vcftools_0.1.13/bin/vcf-merge",
            "bgzip": "/home/Pipeline/tabix/tabix-0.2.6/bgzip",
            "tabix": "/home/Pipeline/tabix/tabix-0.2.6/tabix",
            "#annovar": "",
            "#snpedia": "",
        },
        "databases": {
            "ref": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
            "ref_ucsc_hg19": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
            "ref_mtDNA": "/home/PublicData/h.sapiens_mtDNA/HS_mtDNA.fa",
            "gold_indel": "/home/PublicData/broadinstitute/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
            "oneKG_indel": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf",
            "oneKG_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf",
            "dbsnp": "/home/PublicData/broadinstitute/2.8/hg19/dbsnp_138.hg19.vcf",
            "#dbsnp": "/home/PublicData/dbsnp/human_9606_b150_GRCh37p13/All_20170710.liftover.2.17.11.vcf",
            "onmi_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_omni2.5.hg19.sites.vcf",
            "hapmap_snp": "/home/PublicData/broadinstitute/2.8/hg19/hapmap_3.3.hg19.sites.vcf",
            "target_region": "/home/PublicData/Agilent_v4_71m_reduced.bed",
            "#ensembl_ref_dir": "/home/PublicData/ensembl_GRCh37_75/",
            "#ensembl_ref_fa": "/home/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa",
            "#ensembl_ref_gtf": "/home/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.gtf"
        },
    }
    return default_settings_dict


def mkdir(dir_name):
    os.makedirs(dir_name, exist_ok=True)


def run_pipeline(settings):
    for sample in sorted(settings["samples_dict"]):
        sample_settings = settings
        sample_settings["sample"] = sample
        sample_settings = get_settings_for_SSAP(sample_settings)
        cmd_list = get_cmd_list_for_SSAP(sample_settings)
        write_cmd_list_to_file(sample_settings, cmd_list)
    # debug print
    print(f"# ls {settings['project_script_dir']}")
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do echo $i; done" )
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do qsub $i; done" )
    # ms merge
    run_ms_merge(settings)
    # debug print
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ms.sh ); do qsub $i; done" )


def run_ms_merge(settings):
    project_root = settings["project_root"]
    project_name = os.path.basename(project_root)
    num = len(settings["samples_dict"])
    sample = f"{project_name}_{num}"
    sample_dir = os.path.join(project_root, sample)
    samples_string = " ".join(
        os.path.join(project_root, s, s + ".FINAL.sorted.vcf.gz")
        for s in sorted(settings["samples_dict"])
    )
    ms_WES_merged_vcf = os.path.join(sample_dir, sample + ".merged.vcf")
    cmd_list = [
        f"mkdir -p {sample_dir}",
        f"vcf-merge {samples_string} > {ms_WES_merged_vcf}",
    ]
    script_file = os.path.join(settings["project_script_dir"], sample + ".ms.sh")
    with open(script_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for cmd in cmd_list:
            new_line = cmd + "\n\n"
            f.write(new_line)


def get_settings_for_SSAP(sample_dict):
    sample = sample_dict["sample"]
    sample_dir = os.path.join(sample_dict["project_root"], sample)
    sample_path_prefix = os.path.join(sample_dir, sample)
    _dict = {
        "sample": sample,
        "sample_dir": sample_dir,
        "project_script_dir" : sample_dict["project_script_dir"],
        "number_of_threads": sample_dict["number_of_threads"],

        "bwa": sample_dict["tools"]["bwa"],
        "samtools": sample_dict["tools"]["samtools"],
        "picard": sample_dict["tools"]["picard"],
        "gatk": sample_dict["tools"]["gatk"],
        "vcf_concat": sample_dict["tools"]["vcf_concat"],
        "vcf_sort": sample_dict["tools"]["vcf_sort"],

        "ref": sample_dict["databases"]["ref_ucsc_hg19"],
        "gold_indel": sample_dict["databases"]["gold_indel"],
        "oneKG_indel": sample_dict["databases"]["oneKG_indel"],
        "dbsnp": sample_dict["databases"]["dbsnp"],
        "oneKG_snp": sample_dict["databases"]["oneKG_snp"],
        "hapmap_snp": sample_dict["databases"]["hapmap_snp"],
        "onmi_snp": sample_dict["databases"]["onmi_snp"],

        "read1": sample_dict["samples_dict"][sample]["read1"],
        "read2": sample_dict["samples_dict"][sample]["read2"],

        "sam": sample_path_prefix + ".mem.sam",
        "bam": sample_path_prefix + ".view.bam",
        "sorted_bam": sample_path_prefix + ".sorted.bam",
        "sorted_tmp": sample_path_prefix + ".sorted.tmp",

        "RGID": sample if sample_dict["read_group"]["RGID"] == "__sample__" else sample_dict["read_group"]["RGID"],
        "RGLB": sample if sample_dict["read_group"]["RGLB"] == "__sample__" else sample_dict["read_group"]["RGLB"],
        "RGPL": sample if sample_dict["read_group"]["RGPL"] == "__sample__" else sample_dict["read_group"]["RGPL"],
        "RGPU": sample if sample_dict["read_group"]["RGPU"] == "__sample__" else sample_dict["read_group"]["RGPU"],
        "RGSM": sample if sample_dict["read_group"]["RGSM"] == "__sample__" else sample_dict["read_group"]["RGSM"],
        "RGCN": sample if sample_dict["read_group"]["RGCN"] == "__sample__" else sample_dict["read_group"]["RGCN"],

        "ARRG_bam": sample_path_prefix + ".ARRG.bam",
        "MD_bam": sample_path_prefix + ".MD.bam",
        "MD_metrics_txt": sample_path_prefix + ".MD_picard_metrics.txt",
        "MD_flagstat_txt": sample_path_prefix + ".MD_samtools_flagstat.txt",

        "RTC_intervals": sample_path_prefix + ".RTC.intervals",
        "IR_bam": sample_path_prefix + ".IR.bam",
        "BR_table": sample_path_prefix + ".BR_table.txt",
        "BQSR_BR_table": sample_path_prefix + ".BQSR_BR_table.txt",
        "AC_plot_pdf": sample_path_prefix + ".AC_plot.pdf",
        "BQSR_BR_bam": sample_path_prefix + ".BQSR_BR.bam",

        "gatk_HC_vcf": sample_path_prefix + ".HC.vcf",

        "gatk_VR_SNP_recal": sample_path_prefix + ".HC.VR_SNP.recal",
        "gatk_VR_SNP_tranches": sample_path_prefix + ".HC.VR_SNP.tranches",
        "gatk_VR_SNP_plots": sample_path_prefix + ".HC.VR_SNP_plots.R",
        "gatk_AR_SNP_vcf": sample_path_prefix + ".HC.VQSR_AR_SNP.vcf",

        "gatk_VR_INDEL_recal": sample_path_prefix + ".HC.VR_INDEL.recal",
        "gatk_VR_INDEL_tranches": sample_path_prefix + ".HC.VR_INDEL.tranches",
        "gatk_VR_INDEL_plots": sample_path_prefix + ".HC.VR_INDEL_plots.R",
        "gatk_AR_INDEL_vcf": sample_path_prefix + ".HC.VQSR_AR_SNP.VQSR_AR_INDEL.vcf",

        "gatk_SV_SNP_raw_vcf": sample_path_prefix + ".HC.VQSR.raw_snp.vcf",
        "gatk_SV_SNP_fil_vcf": sample_path_prefix + ".HC.VQSR.fil_snp.vcf",

        "gatk_SV_INDEL_raw_vcf": sample_path_prefix + ".HC.VQSR.raw_indel.vcf",
        "gatk_SV_INDEL_fil_vcf": sample_path_prefix + ".HC.VQSR.fil_indel.vcf",

        "vcftools_concat_vcf": sample_path_prefix + ".FINAL.concat.vcf",
        "vcftools_sorted_vcf": sample_path_prefix + ".FINAL.sorted.vcf",
        "vcftools_sorted_vcf_gz": sample_path_prefix + ".FINAL.sorted.vcf.gz",

        "run_annovar": sample_dict["run_annovar"],
        "annovar_output": sample_path_prefix + ".FINAL.annovar",
        "annovar_txt": sample_path_prefix + ".FINAL.annovar.hg19_multianno.txt",
        "annovar_vcf": sample_path_prefix + ".FINAL.annovar.hg19_multianno.vcf",
        "annovar_header_txt": sample_path_prefix + ".FINAL.annovar.hg19_multianno.header.txt",

        "add_tokens": sample_dict["add_tokens"],
        "debug": sample_dict["debug"],
    }
    return _dict

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

        get_cmd_gatk_RTC_intervals(sample_settings),
        get_cmd_gatk_IR_bam(sample_settings),
        get_cmd_gatk_BR_table(sample_settings),
        get_cmd_gatk_BQSR_BR_table(sample_settings),
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

        get_cmd_bgzip_vcf(sample_settings),
        get_cmd_tabix_gz(sample_settings),

        clear_after_competion(sample_settings)
    ]
    if sample_settings["run_annovar"]:
        cmd_list += [
            get_cmd_annovar(sample_settings),
            get_cmd_annovar_add_header(sample_settings),
        ]
    return cmd_list


###############################################################################
def get_mkdir_cmd(d):
    return "mkdir -p {sample_dir}".format(**d)


###############################################################################
def reduce_spaces_and_newlines(s):
    s = s.replace("\n", " ")
    s = " ".join([i for i in s.split(" ") if i])
    return s


def get_cmd(d):
    d["token"] = "{sample_dir}/token.{sample}.{token_suffix}".format(**d)
    d["flags"] = " && ".join([" [ -f {:s} ] ".format(i) for i in d["files_list"]]) + " && [ ! -f {token} ] ".format(**d)
    if d["add_tokens"]:
        cmd = """
            {flags} &&
            dt1=`date +%y%m%d_%H%M%S` && echo $dt1 {token} &&
            {main_cmd} &&
            du {out_file} > {out_file}.$dt1.du &&
            md5sum {out_file} > {out_file}.$dt1.md5 &&
            dt2=`date +%y%m%d_%H%M%S` &&
            echo $dt1 $dt2 > {token} ||
            echo "TOKEN SKIPPED {token}"
            """.format(**d)
    else:
        cmd = "{main_cmd}".format(**d)
    return reduce_spaces_and_newlines(cmd)


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
    return reduce_spaces_and_newlines(cmd)


###############################################################################
def get_cmd_bwa_mem_sam(d):
    d["token_suffix"] = "bwa_read1_read2_2_sam"
    d["files_list"] = [d["read1"], d["read2"]]
    d["out_file"] = d["sam"]
    d["main_cmd"] = bash_bwa_mem_sam(d)
    return get_cmd(d)


def bash_bwa_mem_sam(d):
    return """{bwa} mem -M -t {number_of_threads} {ref} {read1} {read2} > {sam}""".format(**d)


###############################################################################
def get_cmd_samtools_view_bam(d):
    d["files_list"] = [d["sam"]]
    d["token_suffix"] = "samtools_sam2bam"
    d["out_file"] = d["bam"]
    d["main_cmd"] = bash_samtools_view_bam(d)
    return get_cmd(d)


def bash_samtools_view_bam(d):
    return """{samtools} view -bT {ref} {sam} > {bam}""".format(**d)


###############################################################################
def get_cmd_samtools_sort_bam(d):
    d["files_list"] = [d["bam"]]
    d["token_suffix"] = "samtools_bam_2_sorted_bam"
    d["out_file"] = d["sorted_bam"]
    d["main_cmd"] = bash_samtools_sort_bam(d)
    return get_cmd(d)


def bash_samtools_sort_bam(d):
    return """{samtools} sort -l 9 -O bam -T {sorted_tmp} {bam} > {sorted_bam}""".format(**d)


###############################################################################
def get_cmd_picard_ARRG_bam(d):
    d["files_list"] = [d["sorted_bam"]]
    d["token_suffix"] = "picard_sorted_bam_2_ARRG_bam"
    d["out_file"] = d["ARRG_bam"]
    d["main_cmd"] = bash_picard_ARRG_bam(d)
    return get_cmd(d)


def bash_picard_ARRG_bam(d):
    return """
        {picard}
        AddOrReplaceReadGroups
        INPUT={sorted_bam}
        OUTPUT={ARRG_bam}
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
    d["files_list"] = [d["ARRG_bam"]]
    d["token_suffix"] = "picard_ARRG_bam_2_MD_bam"
    d["out_file"] = d["MD_bam"]
    d["main_cmd"] = bash_picard_MD_bam(d)
    return get_cmd(d)


def bash_picard_MD_bam(d):
    return """
        {picard}
        MarkDuplicates
        INPUT={ARRG_bam}
        OUTPUT={MD_bam}
        METRICS_FILE={MD_metrics_txt}
        ASSUME_SORTED=true
        CREATE_INDEX=true
        VALIDATION_STRINGENCY=LENIENT
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
        """.format(**d)
        # TMP_DIR={tmp_dir}


###############################################################################
def get_cmd_gatk_RTC_intervals(d):
    d["files_list"] = [d["MD_bam"]]
    d["token_suffix"] = "gatk_MD_bam_2_RTC_intervals"
    d["out_file"] = d["RTC_intervals"]
    d["main_cmd"] = bash_gatk_RTC_intervals(d)
    return get_cmd(d)


def bash_gatk_RTC_intervals(d):
    return """
        {gatk}
        -T RealignerTargetCreator
        -R {ref}
        -I {MD_bam}
        -L {target_region}
        -known {gold_indel}
        -known {oneKG_indel}
        -o {RTC_intervals}
        """.format(**d)



###############################################################################
def get_cmd_gatk_IR_bam(d):
    d["files_list"] = [d["MD_bam"], d["RTC_intervals"]]
    d["token_suffix"] = "gatk_MD_bam_RTC_intervals_WE_2_IR_bam"
    d["out_file"] = d["IR_bam"]
    d["main_cmd"] = bash_gatk_IR_bam(d)
    return get_cmd(d)


def bash_gatk_IR_bam(d):
    return """
        {gatk}
        -T IndelRealigner
        -R {ref}
        -I {MD_bam}
        -L {target_region}
        -targetIntervals {RTC_intervals}
        -known {gold_indel}
        -known {oneKG_indel}
        -o {IR_bam}
        """.format(**d)



###############################################################################
def get_cmd_gatk_BR_table(d):
    d["files_list"] = [d["IR_bam"]]
    d["token_suffix"] = "gatk_IR_bam_2_BR_table"
    d["out_file"] = d["BR_table"]
    d["main_cmd"] = bash_gatk_BR_table(d)
    return get_cmd(d)


def bash_gatk_BR_table(d):
    return """
        {gatk}
        -T BaseRecalibrator
        -R {ref}
        -I {IR_bam}
        -knownSites {dbsnp}
        -knownSites {gold_indel}
        -knownSites {oneKG_indel}
        -o {BR_table}
        """.format(**d)


###############################################################################
def get_cmd_gatk_BQSR_BR_table(d):
    d["files_list"] = [d["IR_bam"], d["BR_table"]]
    d["token_suffix"] = "gatk_IR_bam_2_BR_table_BQSR"
    d["out_file"] = d["BQSR_BR_table"]
    d["main_cmd"] = bash_gatk_BQSR_BR_table(d)
    return get_cmd(d)


def bash_gatk_BQSR_BR_table(d):
    return """
        {gatk}
        -T BaseRecalibrator
        -R {ref}
        -I {IR_bam}
        -knownSites {dbsnp}
        -knownSites {gold_indel}
        -knownSites {oneKG_indel}
        -BQSR {BR_table}
        -o {BQSR_BR_table}
        """.format(**d)


###############################################################################
def get_cmd_gatk_AC_pdf(d):
    d["files_list"] = [d["BR_table"], d["BQSR_BR_table"]]
    d["token_suffix"] = "gatk_BR_table_BQSR_BR_table_2_AC_plot_pdf"
    d["out_file"] = d["AC_plot_pdf"]
    d["main_cmd"] = bash_gatk_AC_pdf(d)
    return get_cmd(d)


def bash_gatk_AC_pdf(d):
    return """
        {gatk}
        -T AnalyzeCovariates
        -R {ref}
        -before {BR_table}
        -after {BQSR_BR_table}
        -plots {AC_plot_pdf}
        """.format(**d)


###############################################################################
def get_cmd_gatk_PR_bam(d):
    d["files_list"] = [d["IR_bam"], d["BQSR_BR_table"]]
    d["token_suffix"] = "gatk_IR_bam_BQSR_BR_table_2_BQSR_BR_bam"
    d["out_file"] = d["BQSR_BR_bam"]
    d["main_cmd"] = bash_gatk_PR_bam(d)
    return get_cmd(d)


def bash_gatk_PR_bam(d):
    return """
        {gatk}
        -T PrintReads
        -R {ref}
        -I {IR_bam}
        -BQSR {BQSR_BR_table}
        -o {BQSR_BR_bam}
        """.format(**d)


###############################################################################
def get_cmd_gatk_HC_vcf(d):
    d["files_list"] = [d["BQSR_BR_bam"]]
    d["token_suffix"] = "gatk_BQSR_BR_bam_2_gatk_HC_vcf"
    d["out_file"] = d["gatk_HC_vcf"]
    d["main_cmd"] = bash_gatk_HC_vcf(d)
    return get_cmd(d)


def bash_gatk_HC_vcf(d):
    return """
        {gatk}
        -T HaplotypeCaller
        -R {ref}
        -I {BQSR_BR_bam}
        -D {dbsnp}
        --genotyping_mode DISCOVERY
        -stand_call_conf 30
        -o {gatk_HC_vcf}
        -nct {number_of_threads}
        """.format(**d)


###############################################################################
def get_cmd_gatk_UG_vcf(d): # not ready tested
    d["files_list"] = [d["BQSR_BR_bam"]]
    d["token_suffix"] = "gatk_BQSR_BR_bam_2_gatk_UG_vcf"
    d["out_file"] = d["gatk_UG_vcf"]
    d["main_cmd"] = bash_gatk_UG_vcf(d)
    return get_cmd(d)


def bash_gatk_UG_vcf(d): # not ready tested
    return """
        {gatk}
        -T UnifiedGenotyper
        -R {ref}
        -I {BQSR_BR_bam}
        -D {dbsnp}
        --genotyping_mode DISCOVERY
        -o {gatk_UG_vcf}
        --output_mode EMIT_ALL_SITES
        -nct {number_of_threads}
        --max_num_PL_values 100
        --min_base_quality_score 17
        """.format(**d)


###############################################################################
def get_cmd_gatk_VR_SNP_recal_tranches(d):
    d["files_list"] = [d["gatk_HC_vcf"]]
    d["token_suffix"] = "gatk_HC_vcf_2_VR_SNP_recal"
    d["out_file"] = d["gatk_VR_SNP_recal"]
    d["main_cmd"] = bash_gatk_VR_SNP_recal_tranches(d)
    return get_cmd(d)


def bash_gatk_VR_SNP_recal_tranches(d):
    return """
        {gatk}
        -T VariantRecalibrator
        -R {ref}
        -input {gatk_HC_vcf}
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
        -recalFile {gatk_VR_SNP_recal}
        -tranchesFile {gatk_VR_SNP_tranches}
        -rscriptFile {gatk_VR_SNP_plots}
        """.format(**d)


###############################################################################
def get_cmd_gatk_AR_SNP_vcf(d):
    d["files_list"] = [d["gatk_HC_vcf"], d["gatk_VR_SNP_recal"], d["gatk_VR_SNP_tranches"]]
    d["token_suffix"] = "gatk_HC_vcf_2_AR_SNP_vcf"
    d["out_file"] = d["gatk_AR_SNP_vcf"]
    d["main_cmd"] = bash_gatk_AR_SNP_vcf(d)
    return get_cmd(d)


def bash_gatk_AR_SNP_vcf(d):
    return """
        {gatk}
        -T ApplyRecalibration
        -R {ref}
        -input {gatk_HC_vcf}
        -mode SNP
        --ts_filter_level 99.0
        -recalFile {gatk_VR_SNP_recal}
        -tranchesFile {gatk_VR_SNP_tranches}
        -o {gatk_AR_SNP_vcf}
        """.format(**d)


###############################################################################
def get_cmd_gatk_VR_INDEL_recal_tranches(d):
    d["files_list"] = [d["gatk_AR_SNP_vcf"]]
    d["token_suffix"] = "gatk_AR_SNP_vcf_2_VR_INDEL_recal"
    d["out_file"] = d["gatk_VR_INDEL_recal"]
    d["main_cmd"] = bash_gatk_VR_INDEL_recal_tranches(d)
    return get_cmd(d)


def bash_gatk_VR_INDEL_recal_tranches(d):
    return """
        {gatk}
        -T VariantRecalibrator
        -R {ref}
        -input {gatk_AR_SNP_vcf}
        -resource:mills,known=true,training=true,truth=true,prior=12.0 {gold_indel}
        -an DP
        -an FS
        -an MQRankSum
        -an ReadPosRankSum
        -mode INDEL
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
        --maxGaussians 4
        -recalFile {gatk_VR_INDEL_recal}
        -tranchesFile {gatk_VR_INDEL_tranches}
        -rscriptFile {gatk_VR_INDEL_plots}
        """.format(**d)


###############################################################################
def get_cmd_gatk_AR_INDEL_vcf(d):
    d["files_list"] = [d["gatk_AR_SNP_vcf"], d["gatk_VR_INDEL_recal"], d["gatk_VR_INDEL_tranches"]]
    d["token_suffix"] = "gatk_AR_SNP_vcf_2_AR_SNP_AR_INDEL_vcf"
    d["out_file"] = d["gatk_AR_INDEL_vcf"]
    d["main_cmd"] = bash_gatk_AR_INDEL_vcf(d)
    return get_cmd(d)


def bash_gatk_AR_INDEL_vcf(d):
    return """
        {gatk}
        -T ApplyRecalibration
        -R {ref}
        -input {gatk_AR_SNP_vcf}
        -mode INDEL
        --ts_filter_level 99.0
        -recalFile {gatk_VR_INDEL_recal}
        -tranchesFile {gatk_VR_INDEL_tranches}
        -o {gatk_AR_INDEL_vcf}
        """.format(**d)


###############################################################################
def get_cmd_gatk_SV_SNP_raw_vcf(d):
    d["files_list"] = [d["gatk_AR_INDEL_vcf"]]
    d["token_suffix"] = "gatk_VQSR_vcf_2_SV_SNP_raw_vcf"
    d["out_file"] = d["gatk_SV_SNP_raw_vcf"]
    d["main_cmd"] = bash_gatk_SV_SNP_raw_vcf(d)
    return get_cmd(d)


def bash_gatk_SV_SNP_raw_vcf(d):
    return """
        {gatk}
        -T SelectVariants
        -R {ref}
        -V {gatk_AR_INDEL_vcf}
        -selectType SNP
        -o {gatk_SV_SNP_raw_vcf}
        """.format(**d)


###############################################################################
def get_cmd_gatk_VF_SNP_fil_vcf(d):
    d["files_list"] = [d["gatk_SV_SNP_raw_vcf"]]
    d["token_suffix"] = "gatk_SV_SNP_raw_vcf_2_VF_SNP_fil_vcf"
    d["out_file"] = d["gatk_SV_SNP_fil_vcf"]
    d["main_cmd"] = bash_gatk_VF_SNP_fil_vcf(d)
    return get_cmd(d)


def bash_gatk_VF_SNP_fil_vcf(d):
    return """
        {gatk}
        -T VariantFiltration
        -R {ref}
        -V {gatk_SV_SNP_raw_vcf}
        --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        --filterName "SNP_FAIL"
        -o {gatk_SV_SNP_fil_vcf}
        """.format(**d)


###############################################################################
def get_cmd_gatk_SV_INDEL_raw_vcf(d):
    d["files_list"] = [d["gatk_AR_INDEL_vcf"]]
    d["token_suffix"] = "gatk_VQSR_vcf_2_SV_INDEL_raw_vcf"
    d["out_file"] = d["gatk_SV_INDEL_raw_vcf"]
    d["main_cmd"] = bash_gatk_SV_INDEL_raw_vcf(d)
    return get_cmd(d)


def bash_gatk_SV_INDEL_raw_vcf(d):
    return """
        {gatk}
        -T SelectVariants
        -R {ref}
        -V {gatk_AR_INDEL_vcf}
        -selectType INDEL
        -o {gatk_SV_INDEL_raw_vcf}
        """.format(**d)


###############################################################################
def get_cmd_gatk_VF_INDEL_fil_vcf(d):
    d["files_list"] = [d["gatk_SV_INDEL_raw_vcf"]]
    d["token_suffix"] = "gatk_SV_INDEL_raw_vcf_2_VF_INDEL_fil_vcf"
    d["out_file"] = d["gatk_SV_INDEL_fil_vcf"]
    d["main_cmd"] = bash_gatk_VF_INDEL_fil_vcf(d)
    return get_cmd(d)


def bash_gatk_VF_INDEL_fil_vcf(d):
    return """
        {gatk}
        -T VariantFiltration
        -R {ref}
        -V {gatk_SV_INDEL_raw_vcf}
        --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        --filterName "INDEL_FAIL"
        -o {gatk_SV_INDEL_fil_vcf}
        """.format(**d)


###############################################################################
def get_cmd_vcf_concat(d):
    d["files_list"] = [d["gatk_SV_SNP_fil_vcf"], d["gatk_SV_INDEL_fil_vcf"]]
    d["token_suffix"] = "vcftools_concat_vcf"
    d["out_file"] = d["vcftools_concat_vcf"]
    d["main_cmd"] = bash_vcf_concat(d)
    return get_cmd(d)


def bash_vcf_concat(d):
    return """{vcf_concat} {gatk_SV_SNP_fil_vcf} {gatk_SV_INDEL_fil_vcf} > {vcftools_concat_vcf}""".format(**d)


###############################################################################
def get_cmd_vcf_sort(d):
    d["files_list"] = [d["vcftools_concat_vcf"]]
    d["token_suffix"] = "vcftools_sort_vcf"
    d["out_file"] = d["vcftools_sorted_vcf"]
    d["main_cmd"] = bash_vcf_sort(d)
    return get_cmd(d)


def bash_vcf_sort(d):
    return """cat {vcftools_concat_vcf} | {vcf_sort} > {vcftools_sorted_vcf}""".format(**d)


###############################################################################
def get_cmd_bgzip_vcf(d):
    d["files_list"] = [d["vcftools_sorted_vcf"]]
    d["token_suffix"] = "get_cmd_bgzip_vcf"
    d["out_file"] = d["vcftools_sorted_vcf_gz"]
    d["main_cmd"] = bash_bgzip_vcf(d)
    return get_cmd(d)


def bash_bgzip_vcf(d):
    return """ {bgzip} -c {vcftools_sorted_vcf} > {vcftools_sorted_vcf_gz} """.format(**d)


###############################################################################
def get_cmd_tabix_gz(d):
    d["files_list"] = [d["vcftools_sorted_vcf_gz"]]
    d["token_suffix"] = "get_cmd_tabix_gz"
    d["out_file"] = d["vcftools_sorted_vcf_gz"] + ".tbi"
    d["main_cmd"] = bash_tabix_gz(d)
    return get_cmd(d)


def bash_tabix_gz(d):
    return """ {tabix} {vcftools_sorted_vcf_gz} """.format(**d)


###############################################################################
def get_cmd_annovar(d):
    d["files_list"] = [d["vcftools_sorted_vcf"]]
    d["token_suffix"] = "annovar_vct_and_txt"
    d["out_file"] = d["annovar_txt"]
    d["main_cmd"] = bash_annovar(d)
    return get_cmd(d)

def bash_annovar(d):
    input_file = d["vcftools_sorted_vcf"]
    output_file = d["annovar_output"]
    table_annovar="perl /home/PublicData/annovar_src/annovar_20190101/table_annovar.pl"
    convert2annovar="perl /home/PublicData/annovar_src/annovar_20190101/convert2annovar.pl"
    annovar_db_folder="/home/PublicData/annovar_src/annovar_20190101/humandb"

    db_gene = [
        "refGene",
        "knownGene",
        "ensGene",
    ]

    db_filter = [
        "snp138",
        "avsnp138",
        "avsnp150",
        "ALL.sites.2015_08",
        "AFR.sites.2015_08",
        "AMR.sites.2015_08",
        "SAS.sites.2015_08",
        "EUR.sites.2015_08",
        "EAS.sites.2015_08",
        "esp6500siv2_ea",
        "esp6500siv2_all",
        "esp6500siv2_aa",
        "popfreq_all_20150413",
        "abraom",
        "hrcr1",
        "kaviar_20150923",
        "cg69",
        "dbnsfp35a",
        "dbscsnv11",
        "kgXref",
        "exac03nonpsych",
        "exac03nontcga",
        "gnomad_exome",
        "gnomad_genome",
        "gme",
        "mcap",
        "revel",
        "nci60",
        "icgc21",
        "cosmic68",
        "cosmic70",
        "clinvar_20180603",
        "mitimpact24",
        "regsnpintron",
        "gerp++elem",
        "gerp++gt2",
        "cadd13",
        "fathmm",
        "gwava",
        "eigen",
    ]

    protocol = ",".join(db_gene + db_filter)
    operation = ",".join(["g"]*len(db_gene) + ["f"]*len(db_filter))

    cmd = f"""
        {table_annovar}
        {input_file}
        {annovar_db_folder}
        --protocol {protocol}
        --operation {operation}
        --outfile {output_file}
        --buildver hg19
        --remove
        --otherinfo
        --onetranscript
        --nastring "."
        --vcfinput
    """
    return reduce_spaces_and_newlines(cmd)


###############################################################################
def get_cmd_annovar_add_header(d):
    d["files_list"] = [d["annovar_txt"], d["annovar_vcf"], ]
    d["token_suffix"] = "annovar_add_header"
    d["out_file"] = d["annovar_header_txt"]
    d["main_cmd"] = bash_annovar_add_header(d)
    return get_cmd(d)

def bash_annovar_add_header(d):
    return """ python /home/PublicData/annovar_src/python/add_header.py -v {annovar_vcf} -t {annovar_txt} """.format(**d)


###############################################################################
def write_cmd_list_to_file(sample_settings, cmd_list):
    script_file = os.path.join(
        sample_settings["project_script_dir"],
        sample_settings["sample"] + ".ss.sh",
    )
    with open(script_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for cmd in cmd_list:
            new_line = cmd + "\n\n"
            f.write(new_line)


###############################################################################
if __name__ == "__main__":
    main()
