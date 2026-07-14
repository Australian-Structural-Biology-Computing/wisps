#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Australian-Structural-Biology-Computing/wisps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Australian-Structural-Biology-Computing/wisps
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
workflow WF_WISPS {
    take:
    samplesheet // channel: samplesheet read in from --input
    main:
    tool_set = params.tools.toLowerCase().split(",") as Set
    ch_samplesheet              = samplesheet
    ch_multiqc                  = Channel.empty()
    ch_versions                 = Channel.empty()
    multiqc_report = Channel.empty()
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).first()
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ).first()  : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo ).first()    : Channel.empty()
    ch_multiqc_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_boltz2_aff = ("boltz" in tool_set)
        ? Channel.fromPath(params.boltz2_aff, checkIfExists: true).first()
        : Channel.empty()
    ch_boltz2_conf = ("boltz" in tool_set)
        ? Channel.fromPath(params.boltz2_conf, checkIfExists: true).first()
        : Channel.empty()
    ch_boltz2_mols = ("boltz" in tool_set)
        ? Channel.fromPath(params.boltz2_mols, checkIfExists: true).first()
        : Channel.empty()
    ch_colabfold_alphafold2_params = ("colabfold" in tool_set)
        ? Channel.fromPath(params.colabfold_alphafold2_params, checkIfExists: true).first()
        : Channel.empty()
    ch_alphafold3_params = ("alphafold3" in tool_set)
        ? Channel.fromPath(params.alphafold3_params, checkIfExists: true).first()
        : Channel.empty()
    //
    // WORKFLOW: Run pipeline
    //

        WISPS (
            ch_samplesheet,
            params.mode,
            ch_versions,
            ch_boltz2_aff,
            ch_boltz2_conf,
            ch_boltz2_mols,
            Channel.fromPath(params.colabfold_envdb, checkIfExists: true).first(),
            Channel.fromPath(params.colabfold_uniref30, checkIfExists: true).first(),
            ch_colabfold_alphafold2_params,
            params.colabfold_num_recycles,
            ch_alphafold3_params,
            params.tools,
            ch_multiqc_config,
            ch_multiqc_custom_config,
            ch_multiqc_logo,
            ch_multiqc_methods_description,
            params.outdir,
            params.analysis_batch_size,
            params.colabfold_batch_size,
            params.interaction_neighbours,
            params.pool,
            params.pool_size
        )
    ch_versions = ch_versions.mix(WISPS.out.versions)

    emit:
    multiqc_report = WISPS.out.multiqc_report
    versions   = WISPS.out.versions
    collated_versions = WISPS.out.collated_versions
    ipsae_report = WISPS.out.ipsae_report
    msa_a3m  = WISPS.out.msa_a3m.map{['a3m': it[1]]}
    msa_json = WISPS.out.msa_json.map{['json': it[1]]}
    msa_yaml = WISPS.out.msa_yaml.map{['yaml': it[1]]}
    msa_csv  = WISPS.out.msa_csv.map{['csv': it[1]]}
    
    colabfold_predictions = WISPS.out.colabfold_scores
                                .join(WISPS.out.colabfold_pdb)
                                .map { ['score': params.iptm_threshold > 0 ? new groovy.json.JsonSlurper().parseText(it[1].text).with { (iptm && iptm != 0) ? iptm : ptm } : null,
                                        'pdb': it[2],
                                        'confidence': it[1]] }
                                .filter{params.iptm_threshold == 0 || it.score >= params.iptm_threshold}

    boltz_predictions = WISPS.out.boltz_pae
                            .join(WISPS.out.boltz_cif)
                            .join(WISPS.out.boltz_confidence)
                            .map { ['confidence': it[3], 
                                    'cif': it[2], 
                                    'pae': it[1], 
                                    'score': params.iptm_threshold > 0 ? new groovy.json.JsonSlurper().parseText(it[3].text).with { (iptm && iptm != 0) ? iptm : ptm } : null] }
                            .filter{params.iptm_threshold == 0 || it.score >= params.iptm_threshold}

    af3_predictions = WISPS.out.alphafold3_confidence
                        .join(WISPS.out.alphafold3_cif)
                        .map { ['confidence': it[1], 
                              'cif': it[2], 
                              'score': params.iptm_threshold > 0 ? new groovy.json.JsonSlurper().parseText(it[1].text).with { (iptm && iptm != 0) ? iptm : ptm } : null] }
                        .filter{params.iptm_threshold == 0 || it.score >= params.iptm_threshold}
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
    WF_WISPS (
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
        WF_WISPS.out.multiqc_report
    )
    publish:
    colabfold_predictions = WF_WISPS.out.colabfold_predictions
    boltz_predictions = WF_WISPS.out.boltz_predictions
    af3_predictions = WF_WISPS.out.af3_predictions
    multiqc_report = WF_WISPS.out.multiqc_report // channel: /path/to/multiqc_report.html
    //collated_versions = WF_WISPS.out.collated_versions
    ipsae_report = WF_WISPS.out.ipsae_report
    msa_a3m  = WF_WISPS.out.msa_a3m.ifEmpty([])
    msa_json = WF_WISPS.out.msa_json.ifEmpty([])
    msa_yaml = WF_WISPS.out.msa_yaml.ifEmpty([])
    msa_csv  = WF_WISPS.out.msa_csv.ifEmpty([])

}

output {
    colabfold_predictions {
        path { sample -> {
                sample.pdb >> "colabfold_predictions/pdb/"
                sample.confidence >> "colabfold_predictions/confidence/"        
            }
        }
    }
    boltz_predictions {
        path { sample -> {
                sample.cif >> "boltz_predictions/cif/${sample.cif.name}"
                sample.pae >> "boltz_predictions/pae/${sample.pae.name}"
                sample.confidence >> "boltz_predictions/confidence/${sample.confidence.name}"        
            }
        }
    }
    af3_predictions {
        path { sample -> {
                sample.cif >> "alphafold3_predictions/cif/${sample.cif.name}"
                sample.confidence >> "alphafold3_predictions/confidence/${sample.confidence.name}"        
            }
        }
    }
    multiqc_report {
        path 'multiqc_report.html'
    }
    ipsae_report {
        path "ipsae/"
    }
    msa_a3m  {
        path "mmseqs/a3m/"
        enabled params.save_mmseqs_out
    }
    msa_json  {
        path "mmseqs/json/"
        enabled params.save_mmseqs_out
    }
    msa_yaml {
        path "mmseqs/yaml/"
        enabled params.save_mmseqs_out
    }
    msa_csv  {
        path "mmseqs/csv/"
        enabled params.save_mmseqs_out
    }

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
