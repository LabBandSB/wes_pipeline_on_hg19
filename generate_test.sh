mkdir -p test_fastq_gz
touch test_fastq_gz/sample_01.R1.fastq.gz
touch test_fastq_gz/sample_01.R2.fastq.gz


python wes_pipeline_on_hg19.py

python wes_pipeline_on_hg19.py \
        --draft_settings_json ./__draft_settings__.json \
        --project_root ./test_alignment \
        --script_dir_name scripts \
        --fastq_dirs_list ./test_fastq_gz \
        --sample_delimiter . \
        --fastq_extension .fastq.gz \
        --R1_fastq_extension .R1.fastq.gz \
        --R2_fastq_extension .R2.fastq.gz


python wes_pipeline_on_hg19.py -j test_alignment/scripts/project_settings.json
