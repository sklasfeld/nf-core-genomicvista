#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/genomicvista
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/genomicvista
    Website: https://nf-co.re/genomicvista
    Slack  : https://nfcore.slack.com/channels/genomicvista
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2 // Enable Nextflow DSL2 (Domain Specific Language)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMICVISTA  } from './workflows/genomicvista'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_genomicvista_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_genomicvista_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_GENOMICVISTA {

    take:
    expression_matrix  // channel: [ path(expression_matrix) ]
    tissue_attributes  // channel: [ path(tissue_attributes) ]
    subject_attributes // channel: [ path(subject_attributes) ]
    annotation_gtf     // channel: [ path(annotation_gtf) ]

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GENOMICVISTA (
        expression_matrix,
        tissue_attributes,
        subject_attributes,
        annotation_gtf
    )

    emit:
    results = GENOMICVISTA.out.results // channel: analysis results
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.expression_matrix,
        params.tissue_attributes,
        params.subject_attributes,
        params.annotation_gtf,
        params.chunk_size,
        params.skip_expression_rows,
        params.expression_meta_cols,
        params.filter_tissue,
        params.filter_gene,
        params.filter_transcript,
        params.filter_sex,
        params.filter_age
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_GENOMICVISTA (
        PIPELINE_INITIALISATION.out.expression_data,
        PIPELINE_INITIALISATION.out.tissue_attributes,
        PIPELINE_INITIALISATION.out.subject_attributes,
        PIPELINE_INITIALISATION.out.annotation_gtf,
        PIPELINE_INITIALISATION.out.optional_parameters
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
