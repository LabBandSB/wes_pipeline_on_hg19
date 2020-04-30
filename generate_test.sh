#!/bin/bash

# create test_dir and test_files
mkdir -p test_fastq_gz
touch test_fastq_gz/sample_01.R1.fastq.gz
touch test_fastq_gz/sample_01.R2.fastq.gz

pwd=$(pwd)
# generate __draft_settings__.json # in this file you can edit path to databases
python wes_pipeline_on_hg19.py

# generate project_settings.json # in this file you can edit samples list
python wes_pipeline_on_hg19.py \
  --draft_settings_json ./__draft_settings__.json \
  --project_root ${pwd}/test_alignment \
  --script_dir_name scripts \
  --fastq_dirs_list ${pwd}/test_fastq_gz \
  --sample_delimiter . \
  --fastq_extension .fastq.gz \
  --R1_fastq_extension .R1.fastq.gz \
  --R2_fastq_extension .R2.fastq.gz

# generate individual script files for alignment pipeline
python wes_pipeline_on_hg19.py -j ${pwd}/test_alignment/scripts/project_settings.json
