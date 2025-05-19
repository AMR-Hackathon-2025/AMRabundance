/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UNTAR                  } from '../modules/nf-core/untar/main'
include { MINIMAP2_ALIGNMENT     } from '../subworkflows/local/minimap2_alignment/main'
include { HTSBOX_PILEUP          } from '../modules/local/htsbox/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_amrabundance_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AMRABUNDANCE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:
    ch_ref_data      = Channel.empty()
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    if(params.card_dir){

        // *****************************************
        // Get the structures into a channel.
        // If the folder is compressed, decompress
        // *****************************************
        if(params.card_dir.endsWith('.tar.gz')){

            card_dir = Channel.fromPath(params.card_dir)
                                        .map { it -> [[id: it.baseName],it] }

            UNTAR (card_dir)
                .untar
                .map { meta, dir -> [ file(dir).listFiles() ] }
                .flatten()
                .set{ refs_to_be_mapped }
            ch_versions = ch_versions.mix(UNTAR.out.versions)

        }
        // otherwise, directly use the optional_data within the folder
        else {
            refs_to_be_mapped = Channel.fromPath(params.card_dir+"/**")
        }
    }

    //refs_to_be_mapped.view {i -> i}
    
    refs_to_be_mapped
        .map { it -> [ [ id: it.baseName ], it ] }
        .map { id, seq_id -> [ id, seq_id ] }
        .groupTuple(by: 0)
        .set { ch_ref_data }
        
    //ch_ref_data.view {t -> t}
    //
    // MODULE: Run Minimap2
    //
    MINIMAP2_ALIGNMENT (
        ch_ref_data,
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGNMENT.out.versions.first())

    //
    // MODULE: Run htsbox pileup
    //
    // HTSBOX_PILEUP (
    //     MINIMAP2_ALIGNMENT.out.minimap_align,
    //     refs_to_be_mapped
    // )
    // ch_versions = ch_versions.mix(HTSBOX_PILEUP.out.versions.first())

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
