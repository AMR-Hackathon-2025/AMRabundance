process PRESENCE_ABSENCE_PARSE {
    label 'process_single'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb3700531c7ec639f59f084ab64c05e881d654dcf829db163539f2f0b095e09d/data' :
        'community.wave.seqera.io/library/biopython:1.84--3318633dad0031e7' }"

    input:
    path(tsv)

    output:
    path "presence_absence_summary.tsv", emit: tsv               
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def script_version = 'v1.0.0'
    """
    merge_presence_absence.py . presence_absence_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_presence_absence.py: ${script_version}
    END_VERSIONS
    """
}
