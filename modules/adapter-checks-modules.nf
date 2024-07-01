//
// Discover consensus adapter sequences using bbmerge.sh. 
//
process BBMERGE {
    cpus 8
    

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("adapters.fa"), emit: adapters

    script:
    """
    bbmerge.sh \\
        t=${task.cpus} \\
        in=${reads[0]} \\
        in2=${reads[1]} \\
        outa=adapters.fa
    """
}

//
// Generates a fasta file of completely random sequences to be used as test adapters
//
process RANDOM_SEQUENCE {

    input:
    val n_seqs  
    val seq_length 
    val seed

    output:
    path("random.fa"), emit: adapters

    script:
    """
    #!/usr/bin/env python
    import random
    random.seed(seed)

    bases = ['A', 'T', 'C', 'G']

    def randombase():
        return bases[int(random.random() * 4)]

    def randomseq(length):
        return ''.join([ randombase() for _ in range(length) ])

    with open('random.fa', 'w') as f:
        f.write(''.join([ f'>random_{i+1}\n{randomseq($seq_length)}\n' for i in range($n_seqs) ]))
    """
}

//
// Prefixes the name of the sequence with the sample name so that we can merge them easily.
//
process RENAME_ADAPTERS {
    input:
    tuple val(meta), path(adapters)

    output:
    tuple val(meta), path("${meta}_adapters.fa"), emit: adapters

    script:
    """
    sed 's/Read/${meta}_Read/' ${adapters} > "${meta}_adapters.fa"
    """
}

// 
// Read filtering using bbduk
process BBDUK {

    cpus 16
    publishDir "${params.output}/bbduk", mode: 'symlink'

    input:
    tuple val(meta), path(reads)
    path(adapters)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: passed
    tuple val(meta), path("*_matched.fastq.gz"), emit: matched
    tuple val(meta), path("${meta}_log.txt"),    emit: log
    tuple val(meta), path("${meta}_stats.txt"),  emit: stats

    script:
    """
    bbduk.sh \
        in=${reads[0]} \\
        in2=${reads[1]} \\
        out=${meta}_R1.fastq.gz \\
        out2=${meta}_R2.fastq.gz \\
        outm=${meta}_R1_matched.fastq.gz \\
        outm2=${meta}_R2_matched.fastq.gz \\
        threads=${task.cpus/2} \\
        TODO: ${bbduk_parameters} \\
        ktrim=r mink=11 hdist=2 k=21 tbo ref=$adapters stats=${meta}_stats.txt > ${meta}_log.txt
    """
}


process SUBSAMPLE {
    input:
    tuple val(meta), path(fastq)
    val n // number of reads to extract

    output:
    tuple val(meta), path("*.fastq.gz")    , emit: fastq

    script:
    """
    out=\$(basename $fastq .fastq.gz)
    seqkit \\
        head \\
        -n $n \\
        --threads $task.cpus \\
        $fastq \\
        -o \${out}_subsample.fastq.gz
    """
}

process CUTADAPT {
    cpus 8

    input:
    tuple val(meta), path(reads)
    val adapter
    val args

    output:
    tuple val(meta), path("*trimmed.fastq.gz"), emit: trimmed
    tuple val(meta), path("${meta}_log.txt"),   emit: log

    script:
    """
    cutadapt \
        -a ${adapter} \
        -A ${adapter} \
        ${args} \
        -o ${meta}_R1_trimmed.fastq.gz \
        -p ${meta}_R2_trimmed.fastq.gz \
        -j ${task.cpus} \
        ${reads[0]} ${reads[1]} > ${meta}_log.txt
    """
}