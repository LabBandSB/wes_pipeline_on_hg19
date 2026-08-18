"""Generate per-sample HaloPlex pipeline scripts.

Adapted from ../new_version/generate_new_sh.py with three changes:

1. Reference/tool paths are read from settings.json instead of being hardcoded
   in this file. The WES generator baked /home/PublicData/... into the source,
   which meant editing Python to move a reference file.
2. Templates are resolved relative to this script, not the current working
   directory, so it can be run from anywhere.
3. Sample-name splitting is checked: if the configured delimiter collapses all
   fastq files into one sample, that is a configuration error and we say so
   loudly instead of silently producing a single bogus script.

Usage:
    python3 generate_haloplex_sh.py -j settings.json
"""
import argparse
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))

# Only genuine fallbacks live here — anything path-like belongs in settings.json.
DEFAULTS = {
    "RGPL": "ILLUMINA",
    "RGCN": "NLA",
    "cohort_name": "cohort",
    "interval_padding": "50",
    "platform_unit": "HaloPlex",
    "threads": "4",
    "fastqc": "fastqc",
    "bwa": "bwa",
    "samtools": "samtools",
    "bcftools": "bcftools",
    "java": "java",
    "picard": "picard",
    "gatk": "gatk",
    "vcf_concat": "vcf-concat",
    "vcf_sort": "vcf-sort",
    "vcf_merge": "vcf-merge",
    "bgzip": "bgzip",
    "tabix": "tabix",
}

REQUIRED = [
    "ref", "dbsnp", "target_region", "amplicons_region",
    "cohort_dir", "cohort_name",
    "annovar_dir", "annovar_humandb", "annovar_protocol", "annovar_operation",
]


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

    sample_dict = defaultdict(lambda: defaultdict(str))
    n_fastq = 0
    for fastq in get_files_generator(fastq_dirs_list, fastq_extension):
        n_fastq += 1
        if fastq.endswith(R1_fastq_extension):
            sample = os.path.basename(fastq).split(R1_fastq_delimiter)[0]
            sample_dict[sample]["read1"] = fastq
        elif fastq.endswith(R2_fastq_extension):
            sample = os.path.basename(fastq).split(R2_fastq_delimiter)[0]
            sample_dict[sample]["read2"] = fastq

    # Files named Haloplex_18.i7_160606.R1.fastq.gz split correctly on "." but
    # collapse to a single "Haloplex" if the WES default "_" is kept.
    if n_fastq > 4 and len(sample_dict) <= 1:
        sys.exit(
            f"ERROR: {n_fastq} fastq files collapsed into {len(sample_dict)} sample(s).\n"
            f"       R1_fastq_delimiter={R1_fastq_delimiter!r} is wrong for these filenames.\n"
            f"       For Haloplex_18.i7_160606.R1.fastq.gz the delimiter must be '.'"
        )

    sample_dict = {
        key: value
        for key, value in sample_dict.items()
        if key + "_m" not in sample_dict
    }  # to exclude unmerged samples

    incomplete = [s for s, v in sample_dict.items() if not (v["read1"] and v["read2"])]
    if incomplete:
        print(f"WARNING: {len(incomplete)} sample(s) missing R1 or R2: "
              f"{', '.join(sorted(incomplete)[:10])}", file=sys.stderr)

    return sample_dict


def prepare_sh_for_sample(d, sample_dict):
    d = dict(d)
    d.update(sample_dict)
    sample = d["sample"]
    d["alignment_dir"] = os.path.join(d["project_dir"], sample)

    script_file = os.path.join(d["script_dir"], f"{sample}.sh")
    with open(script_file, "w") as f:
        for line in open(os.path.join(HERE, "bash_sample_template.sh")):
            f.write(line.format(**d))
        for line in open(os.path.join(HERE, "bash_pipeline_template.sh")):
            f.write(line)
    return script_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", required=True)
    args = parser.parse_args()

    settings = json.load(open(args.settings_json))
    d = dict(DEFAULTS)
    d.update({k: v for k, v in settings.items() if not k.startswith("_comment")})

    missing = [k for k in REQUIRED if k not in d]
    if missing:
        sys.exit(f"ERROR: settings.json is missing required keys: {', '.join(missing)}")

    for key in REQUIRED:
        if key.startswith("annovar_") and key in ("annovar_protocol", "annovar_operation"):
            continue
        if not os.path.exists(d[key]):
            print(f"WARNING: {key} does not exist: {d[key]}", file=sys.stderr)

    n_prot = len(d["annovar_protocol"].split(","))
    n_oper = len(d["annovar_operation"].split(","))
    if n_prot != n_oper:
        sys.exit(f"ERROR: annovar_protocol has {n_prot} entries but "
                 f"annovar_operation has {n_oper} — they must match one-to-one.")

    os.makedirs(d["script_dir"], exist_ok=True)
    os.makedirs(d["project_dir"], exist_ok=True)
    os.makedirs(d["cohort_dir"], exist_ok=True)

    samples_dict = load_fastq_samples(settings)
    written = []
    for sample in sorted(samples_dict):
        written.append(prepare_sh_for_sample(d, {
            "sample": sample,
            "read1": samples_dict[sample]["read1"],
            "read2": samples_dict[sample]["read2"],
        }))

    print(f"generated {len(written)} sample scripts in {d['script_dir']}")
    print()
    print("Per-sample stage produces one GVCF each; nothing is genotyped yet.")
    print(f"  ls {d['script_dir']}/*.sh | xargs -P 4 -n 1 bash")
    print("Then, once every sample has finished, genotype the cohort ONCE:")
    print(f"  ./run_cohort.sh {os.path.abspath(args.settings_json)}")


if __name__ == "__main__":
    main()
