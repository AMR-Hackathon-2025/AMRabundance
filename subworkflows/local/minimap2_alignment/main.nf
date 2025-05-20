include { MINIMAP2_INDEX      } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN      } from '../../../modules/nf-core/minimap2/align/main'

workflow MINIMAP2_ALIGNMENT {

    take:

    ch_ref // channel: [meta, ref]
    ch_fasta // channel: [meta2, fasta/fastq]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    MINIMAP2_INDEX ( ch_ref )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    MINIMAP2_INDEX.out.index
        .combine(ch_fasta)
        .multiMap { index_id, index, meta, fastq ->
           reads: [meta, fastq]
           index: [index_id, index]
        }
        .set { ch_mapping_input }

    MINIMAP2_ALIGN ( ch_fasta, MINIMAP2_INDEX.out.index, params.bam_format, params.bam_index_extension, params.cigar_paf_format, params.cigar_bam )
    // MINIMAP2_ALIGN ( ch_mapping_input.reads, ch_mapping_input.index, params.bam_format, params.bam_index_extension, params.cigar_paf_format, params.cigar_bam )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    if (params.bam_format) {
        minimap_out = MINIMAP2_ALIGN.out.bam
    } else {
        minimap_out = MINIMAP2_ALIGN.out.paf
    }

    if (params.bam_index_extension) {
        minimap_index = MINIMAP2_ALIGN.out.index
    } else {
        minimap_index = []
    }
    emit:
    // TODO nf-core: edit emitted channels
    minimap_align = minimap_out       // channel: [ val(meta), [ bam ] ]
    minimap_index = minimap_index     // channel: [ val(meta), [ index ] ]
    versions      = ch_versions       // channel: [ versions.yml ]
}

