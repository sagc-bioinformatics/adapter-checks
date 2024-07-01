#!/usr/bin/env nextflow

params.subsample_n = 100000

// Information on these processes can be found in the module file
include { SUBSAMPLE       } from './modules/adapter-checks-modules.nf'
include { RENAME_ADAPTERS } from './modules/adapter-checks-modules.nf'
include { RANDOM_SEQUENCE } from './modules/adapter-checks-modules.nf'
include { BBMERGE         } from './modules/adapter-checks-modules.nf'
include { BBDUK           } from './modules/adapter-checks-modules.nf'

workflow {
    
    // Take untrimmed read pairs as input
    ch_raw = channel.fromFilePairs(params.input_fastq + "/*_R{1,2}*fastq.gz", size: 2)


    ch_raw_subsample = SUBSAMPLE(ch_raw.transpose(), params.subsample_n).groupTuple(size: 2)
    
    ch_raw_subsample | BBMERGE | RENAME_ADAPTERS

    ch_detected_adapters = RENAME_ADAPTERS.out.adapters.map { it[1] }

    ch_premade_adapter = channel.of(file(params.adapters))

    adapters = ch_premade_adapter
        .mix(ch_detected_adapters)
        .collectFile(name: 'adapters.fa', newLine: false)

    // Now run bbduk on everything
    // ch_trim = channel.fromFilePairs(params.trimmed_dir + "/*_R{1,2}*fastq.gz", size: -1)
    
    // BBDUK(ch_trim, adapters.first())

}
