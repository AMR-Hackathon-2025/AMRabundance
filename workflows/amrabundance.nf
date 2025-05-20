/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UNTAR                                 } from '../modules/nf-core/untar/main'
include { SAMTOOLS_FAIDX                        } from '../modules/nf-core/samtools/faidx/main'
include { MINIMAP2_ALIGNMENT                    } from '../subworkflows/local/minimap2_alignment/main'
include { SAMTOOLS_DEPTH                        } from '../modules/nf-core/samtools/depth/main'
include { PRESENCE_ABSENCE                      } from '../modules/local/calc_presence_absence/main'
include { PRESENCE_ABSENCE_PARSE                } from '../modules/local/calc_presence_absence_parse/main'
include { HTSBOX_PILEUP as HTSBOX_PILEUP_VCF    } from '../modules/local/htsbox/main'
include { HTSBOX_PILEUP as HTSBOX_PILEUP_PILEUP } from '../modules/local/htsbox/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_amrabundance_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AMRABUNDANCE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta      // channel: fasta read in from --fasta
    
    main:
    ch_ref_data      = Channel.empty()
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // if(params.card_dir){

    //     // *****************************************
    //     // Get the structures into a channel.
    //     // If the folder is compressed, decompress
    //     // *****************************************
    //     if(params.card_dir.endsWith('.tar.gz')){

    //         card_dir = Channel.fromPath(params.card_dir)
    //                                     .map { it -> [[id: it.baseName],it] }

    //         UNTAR (card_dir)
    //             .untar
    //             .map { meta, dir -> [ file(dir).listFiles() ] }
    //             .flatten()
    //             .set{ refs_to_be_mapped }
    //         ch_versions = ch_versions.mix(UNTAR.out.versions)

    //     }
    //     // otherwise, directly use the optional_data within the folder
    //     else {
    //         refs_to_be_mapped = Channel.fromPath(params.card_dir+"/**")
    //     }
    // }

    // //refs_to_be_mapped.view {i -> i}
    
    // refs_to_be_mapped
    //     .map { it -> [ [ id: it.baseName ], it ] }
    //     .map { id, seq_id -> [ id, seq_id ] }
    //     .groupTuple(by: 0)
    //     .set { ch_ref_data }
        
    //ch_ref_data.view {t -> t}
    
    SAMTOOLS_FAIDX (
        ch_fasta,
        [ [ id:'no_fai' ],[] ],
        true
    )
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )
    //
    // MODULE: Run Minimap2
    //
    MINIMAP2_ALIGNMENT (
        ch_fasta,
        // ch_ref_data,
        ch_samplesheet
    )
    ch_bam_samtools = MINIMAP2_ALIGNMENT.out.minimap_align
    ch_bam_htsbox   = MINIMAP2_ALIGNMENT.out.minimap_align
    ch_versions     = ch_versions.mix(MINIMAP2_ALIGNMENT.out.versions.first())
    //ch_bam_samtools.view {i -> i}
    //
    // MODULE: Run SAMtools depth
    // 
    ch_bam_samtools
        .map { id, bam -> [ bam ] }
        .map { it -> [ [ id: it.baseName[0] ], it ] }
        .set { ch_samtools_depth }
    
    //ch_samtools_depth.view {i -> i}
    
    SAMTOOLS_DEPTH (
        ch_bam_samtools,
        // ch_samtools_depth,
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())
    
    PRESENCE_ABSENCE (
        SAMTOOLS_DEPTH.out.tsv,
        params.fasta_lengths,
        params.depth_threshold
    )
    ch_versions = ch_versions.mix(PRESENCE_ABSENCE.out.versions.first())
    
    /*
        MODULE: Summarise presence/absence outputs
    */
    PRESENCE_ABSENCE_PARSE (
        PRESENCE_ABSENCE.out.tsv.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(PRESENCE_ABSENCE_PARSE.out.versions)
    
    // ch_bam_htsbox
    //     .map { it -> [ [ id: it.baseName.split("_")[0] ], it ] }
    //     //.map { id -> [ [ ref: id.baseName.split("_")[0] ] ]}
    //     .set { ch_htsbox_input }
    
    // ch_htsbox_input.view {i -> i}

//     ch_samtools_depth
//         .map { id, bam -> [tokens: id.tokenize("_") ]}
//         .set { ch_samtools_depth_htsbox }
    
//     map { id, reads ->
//     (sample, replicate, type) = id.tokenize("_")
//     meta = [sample:sample, replicate:replicate, type:type]
//     [meta, reads]
// }
    
//     ch_samtools_depth_htsbox.view {i -> i}

//     ch_ref_data
//         .combine(ch_bam_htsbox)
//         // .multiMap { ref_id, ref, meta, bam ->
//         //    refs: [ref_id, ref]
//         //    bam: [meta, bam]
//         // }
//         .set { ch_htsbox_input }
    
    //ch_htsbox_input.view {i -> i}
    //
    // MODULE: Run htsbox pileup
    //
    HTSBOX_PILEUP_VCF (
        ch_bam_samtools,
        // ch_samtools_depth,
        ch_fasta,
        SAMTOOLS_FAIDX.out.fai
    )
    ch_versions = ch_versions.mix(HTSBOX_PILEUP_VCF.out.versions.first())

    HTSBOX_PILEUP_PILEUP (
        ch_bam_samtools,
        // ch_samtools_depth,
        ch_fasta,
        SAMTOOLS_FAIDX.out.fai
    )
    ch_versions = ch_versions.mix(HTSBOX_PILEUP_VCF.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'amrabundance_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
