import argparse
import json
import os
from collections import defaultdict

d = {
    "RGID": "__sample__",
    "RGLB": "__sample__",
    "RGPL": "ILLUMINA",
    "RGPU": "SureSelectV4",
    "RGSM": "__sample__",
    "RGCN": "NLA",

    "fastqc": "fastqc",
    "bwa": "bwa",
    "samtools": "samtools",
    "bcftools": "bcftools",
    "java": "java",
    "picard": "picard",
    "gatk": "gatk3",
    "vcf_concat": "vcf-concat",
    "vcf_sort": "vcf-sort",
    "vcf_merge": "vcf-merge",
    "bgzip": "bgzip",
    "tabix": "tabix",

    "dbsnp": "/home/PublicData/broadinstitute/2.8/hg19/dbsnp_138.hg19.vcf",
    "gold_indel": "/home/PublicData/broadinstitute/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
    "hapmap_snp": "/home/PublicData/broadinstitute/2.8/hg19/hapmap_3.3.hg19.sites.vcf",
    "oneKG_indel": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf",
    "oneKG_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf",
    "onmi_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_omni2.5.hg19.sites.vcf",
    "ref": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
    "ref_mtDNA": "/home/PublicData/h.sapiens_mtDNA/HS_mtDNA.fa",
    "ref_ucsc_hg19": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
    "target_region": "/home/PublicData/Agilent_v4_71m_reduced.bed",
}


def get_files_generator(dirs_list, extension=""):
    for path in dirs_list:
        for data_file in os.listdir(path):
            if data_file:
                data_path = os.path.join(path, data_file)
                if os.path.isfile(data_path) and data_path.endswith(extension):
                    yield data_path
                elif os.path.isdir(data_path):
                    yield from get_files_generator([data_path], extension)


def load_fastq_samples(settings):
    fastq_dirs_list = settings["fastq_dirs_list"]
    fastq_extension = settings["fastq_extension"]
    R1_fastq_extension = settings["R1_fastq_extension"]
    R2_fastq_extension = settings["R2_fastq_extension"]
    R1_fastq_delimiter = settings["R1_fastq_delimiter"]
    R2_fastq_delimiter = settings["R2_fastq_delimiter"]
    #
    sample_dict = defaultdict(lambda: defaultdict(str))
    for fastq in get_files_generator(fastq_dirs_list, fastq_extension):
        if fastq.endswith(R1_fastq_extension):
            sample = os.path.basename(fastq).split(R1_fastq_delimiter)[0]
            sample_dict[sample]["read1"] = fastq
        elif fastq.endswith(R2_fastq_extension):
            sample = os.path.basename(fastq).split(R2_fastq_delimiter)[0]
            sample_dict[sample]["read2"] = fastq
    sample_dict = {
        key: value
        for key, value in sample_dict.items()
        if key + "_m" not in sample_dict
    }  # to exclude unmerged samples
    #
    return sample_dict


def samples_dict_to_list(d):
    arr = [
        {
            'sample': sample,
            'read1': d[sample]['read1'],
            'read2': d[sample]['read2'],
        } for sample in d
    ]
    return arr


def prepare_sh_for_sample(sample_dict):
    d.update(sample_dict)
    sample = d["sample"]
    project_dir = d["project_dir"]
    script_dir = d['script_dir']

    d["alignment_dir"] = os.path.join(project_dir, sample)

    script_name = f"{sample}.sh"
    script_file = os.path.join(script_dir, script_name)
    with open(script_file, "w") as f:
        for line in open("bash_sample_template.sh"):
            new_line = line.format(**d)
            f.write(new_line)
        for line in open("bash_pipeline_template.sh"):
            f.write(line)


def parse_arguments_to_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", default=None, required=True)
    args = parser.parse_args()
    if args.settings_json:
        settings = json.load(open(args.settings_json))
    else:
        settings = [{}]
    return settings


def main():
    settings = parse_arguments_to_settings()
    samples_dict = load_fastq_samples(settings)
    samples_list = samples_dict_to_list(samples_dict)

    script_dir = settings['script_dir']
    os.makedirs(script_dir, exist_ok=True)

    for sample_dict in samples_list:
        sample_dict['script_dir'] = settings['script_dir']
        sample_dict['project_dir'] = settings['project_dir']
        prepare_sh_for_sample(sample_dict)


if __name__ == '__main__':
    main()
