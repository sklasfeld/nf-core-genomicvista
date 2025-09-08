/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_genomicvista_pipeline'

// Local modules
include { DATA_PREPROCESSING    } from '../modules/local/data_preprocessing'
include { EXPRESSION_STATISTICAL   } from '../modules/local/expression_statistical'
include { VISUALIZATION         } from '../modules/local/visualization'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMICVISTA {

    take:
    ch_expression_matrix     // channel: [ path(expression_matrix) ]
    ch_tissue_attributes   // channel: [ path(tissue_attributes) ]
    ch_subject_attributes  // channel: [ path(subject_attributes) ]
    ch_annotation_gtf     // channel: [ path(annotation_gtf) ]

    main:

    ch_versions = Channel.empty()


    //
    // MODULE: Data preprocessing 
    //
    DATA_PREPROCESSING (
        ch_expression_matrix,
        ch_tissue_attributes,
        ch_subject_attributes,
        ch_annotation_gtf
    )
    ch_versions = ch_versions.mix(DATA_PREPROCESSING.out.versions.first())

    //
    // MODULE: Statistical analysis of expression data
    //
    EXPRESSION_STATISTICAL (
        DATA_PREPROCESSING.out.processed_data
    )
    ch_versions = ch_versions.mix(EXPRESSION_STATISTICAL.out.versions.first())


    //
    // MODULE: Visualization
    //
    VISUALIZATION (
        EXPRESSION_STATISTICAL.out.analysis_results
    )
    ch_versions = ch_versions.mix(VISUALIZATION.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'genomicvista_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    processed_data    = DATA_PREPROCESSING.out.processed_data
    analysis_results  = EXPRESSION_STATISTICAL.out.analysis_results
    plots            = VISUALIZATION.out.plots
    results          = VISUALIZATION.out.plots.collect()
    versions       = ch_collated_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
