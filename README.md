

# adapter-checks

A pipeline to validate read trimming preformance. Assumes paired end.
Dependencies are the latest stable version of `nextflow` and `conda` (alternatively, BYO versions of the tools found in `resources/environment.yml` and remove all `conda` configuration in `nextflow.config`).

1. `cutadapt` is used to trim reads, given a known adapter sequence (`--cutadapt_adapter`). Other arguments can be provided with `--cutadapt_args`.
2. Independently, `bbmerge` is used to automatically detect adapter sequences from untrimmed data. Untrimmed files are subsetted to a fixed number of reads, configurable with `--subsample_n`. 
3. `bbmerge` consensus sequences are compared to the adapter provided to `cutadapt`. (TODO)
4. Consensus sequences are merged with a premade adapter list (can be altered via `--adapters`). The default uses the list provided by `bbmap`.
5. `bbduk` is used with the merged adapter list to filter for adapters. Two different sets of arguments are used, `-bbduk_standard_args` and `-bbduk_lenient_args`.
6. Random adapter sequences are generated and the same parameters are used with the random adapters as a comparison.
7. Report is generated (TODO). Most frequent matches are compared across each parameter combination.

## Running

```{bash}
WORKFLOW=path/to/adapter-checks/

nextflow run path/to/adapter-checks \
    --input_fastq input_dir \
    --cutadapt_adapter 'CTGTCTCTTATACACATCT' \
    --cutadapt_args '-e 0.3' \
    --bbduk_standard_args 'ktrim=r mink=11 hdist=2 k=21 tbo' \
    --bbduk_lenient_args 'ktrim=r mink=11 hdist=1 k=23 tbo' \
    --output results_dir \
    -resume
```

## Output

A report will be generated with some custom figures.

TODO: It would be nice to have those as part of a multiQC plugin, perhaps.