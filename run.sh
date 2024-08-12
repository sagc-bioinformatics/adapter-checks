module load nextflow

WORKFLOW=path/to/adapter-checks/

nextflow run path/to/adapter-checks \
    --merged_fastq merged \
    --adapters $WORKFLOW/resources/bbmap_adapters.fa \
    --bbduk_standard_args 'ktrim=r mink=11 hdist=2 k=21 tbo' \
    --bbduk_lenient_args 'ktrim=r mink=11 hdist=1 k=23 tbo' \
    --output results_dir \
    --conda_env $WORKFLOW/resources/environment.yml \
    -resume

