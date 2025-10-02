#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/wisps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/wisps
    Website: https://nf-co.re/wisps
    Slack  : https://nfcore.slack.com/channels/wisps
----------------------------------------------------------------------------------------
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { WISPS  } from './workflows/wisps'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_wisps_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_wisps_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_WISPS {
    take:
    samplesheet // channel: samplesheet read in from --input
    main:
    ch_samplesheet              = samplesheet
    ch_multiqc                  = Channel.empty()
    ch_versions                 = Channel.empty()
    multiqc_report = Channel.empty()
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).first()
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ).first()  : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo ).first()    : Channel.empty()
    ch_multiqc_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    //
    // WORKFLOW: Run pipeline
    //
    
        WISPS (
            ch_samplesheet,
            params.mode,
            ch_versions,
            Channel.fromPath(params.boltz_ccd, checkIfExists: true).first(),
            Channel.fromPath(params.boltz_model, checkIfExists: true).first(),
            Channel.fromPath(params.boltz2_aff, checkIfExists: true).first(),
            Channel.fromPath(params.boltz2_conf, checkIfExists: true).first(),
            Channel.fromPath(params.boltz2_mols, checkIfExists: true).first(),
            Channel.fromPath(params.colabfold_db, checkIfExists: true).first(),
            Channel.fromPath(params.uniref30, checkIfExists: true).first(),
            params.mmseqs_batch_size,
            params.colabfold_model_preset,
            Channel.fromPath(params.colabfold_alphafold2_params, checkIfExists: true).first(),
            params.colabfold_num_recycles,
            Channel.fromPath(params.alphafold3_params, checkIfExists: true).first(),
            params.tools,
            ch_multiqc_config,
            ch_multiqc_custom_config,
            ch_multiqc_logo,
            ch_multiqc_methods_description,
            params.outdir
        )
        ch_versions = ch_versions.mix(WISPS.out.versions)
  
    emit:
    multiqc_report //= WISPS.out.multiqc_report // channel: /path/to/multiqc_report.html
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
        params.input
    )
    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_WISPS (
        PIPELINE_INITIALISATION.out.samplesheet
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
        NFCORE_WISPS.out.multiqc_report
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
