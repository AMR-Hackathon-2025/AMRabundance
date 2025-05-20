process HTSBOX_PILEUP {
    tag "$meta.id"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b615c176495c90c5ed3c17c8983e71b521cbb61b6ebca5b9457b736ff053370d/data' :
        'community.wave.seqera.io/library/htsbox:r346--203bdf1cfada0cc6' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf")                       , optional: true, emit: vcf
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def htsbox_version = 'r346'
    """
    htsbox pileup \\
        $args \\
        -f $reference \\
        $bam > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htsbox: ${htsbox_version}
    END_VERSIONS
    """
}
