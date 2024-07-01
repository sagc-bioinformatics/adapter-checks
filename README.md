

# adapter-checks

A pipeline to validate read trimming preformance. Assumes paired end.

1. `cutadapt` is used to trim reads, given a known adapter sequence (`--cutadapt_adapter`). Other arguments can be provided with `--cutadapt_args`.
2. Independently, `bbmerge` is used to automatically detect adapter sequences from untrimmed data. Untrimmed files are subsetted to a fixed number of reads, configurable with `--subsample_n`. 
3. `bbmerge` consensus sequences are compared to the adapter provided to `cutadapt`. (TODO)
4. Consensus sequences are merged with a premade adapter list (provided through `--adapters`). The example here uses the list provided by `bbmap`.
5. `bbduk` is used with the merged adapter list to filter for adapters. Two different sets of arguments are used, `-bbduk_standard_args` and `-bbduk_lenient_args`.
6. Report is generated (TODO)


## Running

```{bash}
WORKFLOW=path/to/adapter-checks/

nextflow run path/to/adapter-checks \
    --input_fastq input_dir \
    --cutadapt_adapter 'CTGTCTCTTATACACATCT' \
    --cutadapt_args '-e 0.3' \
    --adapters $WORKFLOW/resources/bbmap_adapters.fa \
    --bbduk_standard_args 'ktrim=r mink=11 hdist=2 k=21 tbo' \
    --bbduk_lenient_args 'ktrim=r mink=11 hdist=1 k=23 tbo' \
    --output results_dir \
    --conda_env $WORKFLOW/resources/environment.yml \
    -resume
```