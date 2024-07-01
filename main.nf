#!/usr/bin/env nextflow

params.subsample_n = 100000

// Information on these processes can be found in the module file
include { SUBSAMPLE       } from './modules/adapter-checks-modules.nf'
include { RENAME_ADAPTERS } from './modules/adapter-checks-modules.nf'
include { RANDOM_SEQUENCE } from './modules/adapter-checks-modules.nf'
include { BBMERGE         } from './modules/adapter-checks-modules.nf'
include { BBDUK           } from './modules/adapter-checks-modules.nf'
include { CUTADAPT        } from './modules/adapter-checks-modules.nf'
include { SEQ_COMPARE     } from './modules/adapter-checks-modules.nf'

workflow {
    
    // Take untrimmed read pairs as input
    ch_raw = channel.fromFilePairs(params.input_fastq + "/*_R{1,2}*fastq.gz", size: 2)

    // Perform adapter trimming using cutadapt, with given adapter sequence
    /*
    CUTADAPT (
        ch_raw,
        params.cutadapt_adapter,
        params.cutadapt_args,
    )*/

    // Subsample the reads before we use bbmerge (purely for efficiency)
    ch_raw_subsample = SUBSAMPLE(ch_raw.transpose(), params.subsample_n).groupTuple(size: 2)
    
    // Detect adapter sequence with bbmerge, and rename the output sequences so they can all go in one file    
    ch_raw_subsample | BBMERGE | RENAME_ADAPTERS

    // Get premade adapter fasta file
    ch_premade_adapter = channel.of(file(params.adapters))
    
    // Get list of detected adapters (one per sample) 
    ch_detected_adapters = RENAME_ADAPTERS.out.adapters.map { it[1] }

    // Check whether provided adapter was detected by bbmerge
    SEQ_COMPARE (
        params.cutadapt_adapter,
        ch_detected_adapters
    )
    // Print out the percentage of detected sequences which start with the sequence given to cutadapt
    SEQ_COMPARE.out.map {
        n_seq, n_match -> {
            Float.valueOf(n_match.text) / Float.valueOf(n_seq.text)
        }
    }
    .mean()
    .map {
        println "Detected adapter sequences which match input adapter: ${it*100}%"
    }

    // Mix all adapter sequences into one file
    ch_combined_adapters = ch_premade_adapter
        .mix(ch_detected_adapters)
        .collectFile(name: 'adapters.fa', newLine: false)
        .first()

    // Get number of sequences in combined adapters file
    ch_n_seq = ch_combined_adapters | splitFasta | collect | map { it.size() }

    // Create an equivalently long fasta file with purely random sequences
    RANDOM_SEQUENCE (
        ch_n_seq, // number of sequences to be generated
        50, // length in bp
        100 // seed for random generation
    )

    /*
    // Use bbduk to filter for adapters with the two different sets
    BBDUK (
        ch_trim,
        ch_combined_adapters.mix(RANDOM_SEQUENCE.out.adapters)
    )
    */
}
