import os

from generate_new_sh import d

bin_targets = [
    "fastqc",
    "bwa",
    "samtools",
    "bcftools",
    "java",
    "picard",
    "gatk",
    "vcf_concat",
    "vcf_sort",
    "vcf_merge",
    "bgzip",
    "tabix",
]

db_targets = [
    "dbsnp",
    "gold_indel",
    "hapmap_snp",
    "oneKG_indel",
    "oneKG_snp",
    "onmi_snp",
    "ref",
    "ref_mtDNA",
    "ref_ucsc_hg19",
    "target_region",
]


def main():
    for bin in bin_targets:
        _ = os.system(f"which {d[bin]}")

    for db in db_targets:
        if os.path.isfile(d[db]):
            print(f"EXISTS {db}={d[db]}")
        else:
            print(f"NOT EXISTS {db}={d[db]}")


if __name__ == "__main__":
    main()
