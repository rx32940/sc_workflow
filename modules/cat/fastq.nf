process CAT_FASTQ {
    tag "$meta.id"
    label 'process_high_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz")

    script:
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        if (readList.size >= 1) {
            """
            cat ${readList.join(' ')} > ${prefix}.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        }
    } else {
        if (readList.size >= 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }

            """
            cat ${read1.join(' ')} > ${prefix}_S1_L001_R1_001.fastq.gz 
            cat ${read2.join(' ')} > ${prefix}_S1_L001_R2_001.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        }
    }
}