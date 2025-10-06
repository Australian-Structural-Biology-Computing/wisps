/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { PREPARE_INTERACTIONS } from '../modules/local/prepare_interactions'
include { COLLECT_CONFIDENCE } from '../modules/local/collect_confidence'
//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_wisps_pipeline'
include { MSA } from '../subworkflows/local/msa'

//
// LOCAL MODULE: Boltz
//
include { RUN_BOLTZ } from '../modules/local/run_boltz'
include { RUN_ALPHAFOLD3 } from '../modules/local/run_alphafold3'
include { COLABFOLD_BATCH } from '../modules/local/colabfold_batch'
include { BOLTZ_FASTA } from '../modules/local/data_convertor/boltz_fasta'
include { SPLIT_MSA } from '../modules/local/msa_manager/split_msa'
include { IPSAE } from '../modules/local/ipsae'


//
// FUNCTIONS
//
include { getFastaSequences   } from '../subworkflows/local/msa'
include { getYamlSequences    } from '../subworkflows/local/msa'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WISPS {
    
    take:
    ch_samplesheet  // channel: samplesheet read from --input
    mode
    ch_versions     // channel: [ path(versions.yml) ]
    ch_boltz_ccd    // channel: [ path(boltz_ccd) ]
    ch_boltz_model  // channel: [ path(model) ]
    ch_boltz2_aff   // channel: [ path(boltz2_aff) ]
    ch_boltz2_conf  // channel: [ path(boltz2_conf) ]
    ch_mols         // channel: [ path(mols) ]
    ch_colabfold_db // channel: [ path(colabfold_db) ]
    ch_uniref30     // channel: [ path(uniref30) ]
    mmseqs_batch_size // number
    colabfold_model_preset
    ch_colabfold_params
    num_recycles
    ch_alphafold3_params // channel: path(alphafold2_params)
    tools
    ch_multiqc_config
    ch_multiqc_custom_config
    ch_multiqc_logo
    ch_multiqc_methods_description
    outdir
    analysis_batch_size

    main:
    ch_multiqc_files = Channel.empty()
    ch_confidence_scores = Channel.empty()
    

    interaction_mode = mode.split(",").collect { pair ->
        pair.split('-').sort().join('-')
    }
    if ("all-all" in interaction_mode) {
        interaction_mode = ["all-all"]
    }else{
        required = interaction_mode.collect{it.split("-")}.flatten().unique()
        ch_samplesheet
            .map{it[0].group}
            .unique()
            .toList()
            .subscribe{
                emitted -> {
                    if (required - emitted){
                        log.error "The groups ${required - emitted} is not in the sample sheet!"
                        exit 1      
                    }
                }
            }
    }
    log.info "Used Interactions (Sorted): ${interaction_mode.join(',')}"
    ch_samplesheet.map{it[0]}
        .combine(ch_samplesheet.map{it[0]})
        .map{[
            "id": [it[0]["id"], it[1]["id"]].min() + "-" + [it[0]["id"], it[1]["id"]].max(),
            "group": [it[0]["group"], it[1]["group"]].min() + "-" + [it[0]["group"], it[1]["group"]].max(),
            ]}
        .filter{ it["group"] in interaction_mode || interaction_mode[0] == "all-all"}
        .unique()
        .set{ch_unique_pairs}

    ch_samplesheet
    .combine(ch_samplesheet)
    .map{[
            [  
                "id": [it[0]["id"], it[2]["id"]].min() + "-" + [it[0]["id"], it[2]["id"]].max(),
                "group": [it[0]["group"], it[2]["group"]].min() + "-" + [it[0]["group"], it[2]["group"]].max(),
            ],
            [it[0], it[2]],
            [it[1], it[3]],
    ]}.set{ch_pairs}
    
    ch_pairs
    .join(ch_unique_pairs)
    .set{ch_interaction_in}

    //ch_interaction_in.view()
    
    ch_interaction_in
    .filter{it[1][0].type == "protein" && it[1][1].type == "protein"}
    .collectFile{
        seqs = []
        for (j in [0, 1]){
            if (it[2][j].extension == "fasta" || it[2][j].extension == "fa")
                {seqs.add(getFastaSequences(it[2][j].text))}
            else if (it[2][j].extension == "yaml" || it[2][j].extension == ".yml")
                {seqs.add(getYamlSequences(it[2][j].text))}
            else if (it[2][j].extension == "json")
                {seqs.add("")}
        }
        [ "${it[0].id}.fasta", ">${it[0].id}\n${seqs.collect{it.sequence.join(":")}.join(":")}\n" ]
    }
    .map{[["id": it.baseName], it]}
    .set {ch_protein_pairs}

    ch_interaction_in
    .map{
        if (it[1][0].type != "protein" && it[1][1].type == "protein"){
            [it[1][1], it[2][1]]   
        }else if (it[1][1].type != "protein" && it[1][0].type == "protein"){
            [it[1][0], it[2][0]]   
        }
    }.unique()
    .set{ch_single_protein}
    
    //ch_protein_pairs.view()
    MSA (
        ch_protein_pairs
        .mix(
            ch_single_protein
        ),
        ch_colabfold_db,
        ch_uniref30,
        mmseqs_batch_size
    )
    ch_versions = ch_versions.mix(MSA.out.versions)
    
    //MSA.out.json.view()
    //ch_interaction_in.view()
    
    // Prepare interactions for boltz
    ch_boltz_data = Channel.empty()
    ch_split_msa_in = Channel.empty()
    ch_boltz_interactions_in = Channel.empty()

    if ("boltz" in tools.split(",")){
        ch_split_msa_in = MSA.out.a3m.join(ch_protein_pairs.map{it[0]})
        ch_boltz_interactions_in = ch_interaction_in  
    }

    // Adding non protein for boltz
    ch_boltz_data = ch_boltz_data.mix(ch_boltz_interactions_in
        .filter{it[1][0].type != "protein" && it[1][1].type != "protein"}
        .map{[it, ["", ""]]}
    )

    // Adding single protein for boltz
    ch_boltz_data = ch_boltz_data.mix(
        ch_boltz_interactions_in
        .filter{[it[1][0].type, it[1][1].type].count("protein") == 1}
        .map{
            if (it[1][0].type == "protein"){
                [it[1][0], it]
            }else if (it[1][1].type == "protein"){
                [it[1][1], it]
            }
        }
        .combine(MSA.out.a3m.filter{it[0].type == "protein"})
        .filter{it[0] == it[2]}
        .map{
            if (it[1][1][0].type == "protein")
                [it[1], [it[3], ""]]
            else
                [it[1], ["", it[3]]]
        }
    )
    split_batch_id = 0
    SPLIT_MSA(
        ch_split_msa_in
        .map{it[1]}
        .buffer( size: mmseqs_batch_size, remainder: true )
        .map{ split_batch_id += 1; [["id" : "batch-${split_batch_id}-${it.size()}"], it]}
    )

    ch_versions = ch_versions.mix(SPLIT_MSA.out.versions)
    
    SPLIT_MSA.out.msa_csv
    .map{it[1]}
    .flatten()
    .map{[["id" : it.baseName.split("_")[0]], it]}
    .groupTuple(sort: true)
    .set{ch_split_msa_out}

    ch_boltz_data = ch_boltz_data
        .mix(ch_boltz_interactions_in
        .filter{it[1][0].type == "protein" && it[1][1].type == "protein"}
        .map{[["id": it[0].id], it]}
        .join(
            ch_split_msa_out
            .map{
                if (it[1].size() == 1){
                    it[1] = [it[1], it[1]] 
                }
                it
            }
            
        )
        .map{[it[1], it[2]]}
    )
 
    batch_id = 0

    PREPARE_INTERACTIONS(
        ch_boltz_data
            .map{it[0][0].id}
            .buffer( size: analysis_batch_size, remainder: true )
            .map{ batch_id += 1; [["id" : "batch-${batch_id}-${it.size()}"], it]},
        ch_boltz_data
            .map{it[0][2].collect{it.name}}
            .buffer( size: analysis_batch_size, remainder: true ),
        ch_boltz_data
            .map{it[1].collect{it ? it.name: ""}}
            .buffer( size: analysis_batch_size, remainder: true ),
        ch_boltz_data
            .map{[it[0][1][0].type, it[0][1][1].type]}
            .buffer( size: analysis_batch_size, remainder: true ),
        ch_boltz_data
            .map{[it[0][2][0], it[0][2][1], it[1][0], it[1][1]]}
            .buffer( size: analysis_batch_size, remainder: true )
            .map{it.flatten()findAll{it}.unique{it.toUriString()}}
        
    )

    PREPARE_INTERACTIONS.out.fasta
    .map{it[1]}
    .flatten()
    .map{[it.baseName.replace("_boltz_interaction_input", ""), it]}
    .join(
        ch_boltz_data
        .map{[it[0][0].id, [it[1][0], it[1][1]]]}
    )
    .buffer( size: analysis_batch_size, remainder: true )
    .set{ch_boltz_in}

    //ch_boltz_in.view()
    boltz_batch = 0
    
    RUN_BOLTZ(
        ch_boltz_in.map{boltz_batch += 1; [["id": "batch-${boltz_batch}-${it.size()}"], it.collect{it[1]}]},
        ch_boltz_in.map{it.collect{it[2]}.flatten().findAll{it}.unique{it.toUriString()}},
        ch_boltz_model,
        ch_boltz_ccd,
        ch_boltz2_aff,
        ch_boltz2_conf,
        ch_mols
    )
    

    RUN_BOLTZ.out.confidence
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[1], "model": "boltz"], it]}
    .set{ch_boltz_confidence}
    
    RUN_BOLTZ.out.pae
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[1], "model": "boltz"], it]}
    .set{ch_boltz_pae}
    
    RUN_BOLTZ.out.cif
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[0], "model": "boltz"], it]}
    .set{ch_boltz_cif}

    ch_boltz_confidence
    .collect(flat: false, sort: true).multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_boltz_confidence_scores}
   
    //prepare interactions for colabfold
    ch_colabfold_interaction_in = Channel.empty()
    if ("colabfold" in tools.split(",")){
        colabfold_batch = 0
        MSA.out.a3m.join(ch_protein_pairs.map{it[0]})
        .map{it[1]}
        .buffer( size: analysis_batch_size, remainder: true )
        .map{
            colabfold_batch += 1;
            [["id": "batch-${colabfold_batch}"], it]
        }
        .set{
            ch_colabfold_interaction_in
        }
    }
    //ch_colabfold_interaction_in.view()
    
    COLABFOLD_BATCH(
        ch_colabfold_interaction_in,
        colabfold_model_preset,
        ch_colabfold_params,
        [],
        [],
        num_recycles
    )
    ch_versions = ch_versions.mix(COLABFOLD_BATCH.out.versions)
   
    COLABFOLD_BATCH.out.top_ranked_scores
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[0], "model": "colabfold"], it]}
    .set{ch_colabfold_scores}
    
    COLABFOLD_BATCH.out.pdb
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[0], "model": "colabfold"], it]}
    .set{ch_colabfold_pdb}


    ch_colabfold_scores
    .collect(flat: false, sort: true).multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_colabfold_confidence_scores}

    
    af3_batch = 0
    ch_alphafold3_interaction_in = Channel.empty()
    if ("alphafold3" in tools.split(",")){
        MSA.out.a3m.join(ch_protein_pairs.map{it[0]})
        .map{it[1]}
        .buffer( size: analysis_batch_size, remainder: true )
        .map{
            af3_batch += 1;
            [["id": "batch-${af3_batch}-${it.size()}"], it]
        }
        .set{
            ch_alphafold3_interaction_in
        }
    }

    RUN_ALPHAFOLD3 (
        ch_alphafold3_interaction_in,
        ch_alphafold3_params
    )
    ch_versions = ch_versions.mix(RUN_ALPHAFOLD3.out.versions)

    RUN_ALPHAFOLD3.out.top_ranked_cif
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[0], "model": "alphafold3"], it]}
    .set{ch_alphafold3_cif}
    
    RUN_ALPHAFOLD3.out.confidence
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[0], "model": "alphafold3"], it]}
    .set{ch_alphafold3_confidence}
    
    RUN_ALPHAFOLD3.out.summary_confidences
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_")[0], "model": "alphafold3"], it]}
    .set{ch_alphafold3_summary_confidences}


    ch_alphafold3_summary_confidences.collect(flat: false, sort: true).multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_alphafold3_confidence_scores}

    
    ch_colabfold_scores
    .join(ch_colabfold_pdb)
    .map{[["id": it[0].id], it, []]}
    .join(ch_interaction_in
        .filter{it[1][0].type == "protein" || it[1][1].type == "protein"}
        .map{["id": it[0].id]}
    )
    .mix(
        ch_boltz_pae
        .join(ch_boltz_cif)
        .join(ch_boltz_confidence)
        .map{[["id": it[0].id], [it[0], it[1], it[2]], it[3]]}
        .join(ch_interaction_in
            .filter{it[1][0].type == "protein" || it[1][1].type == "protein"}
            .map{["id": it[0].id]}
        )
    )
    .mix(
        ch_alphafold3_confidence.join(ch_alphafold3_cif)
        .map{[["id": it[0].id], it, []]}
        .join(ch_interaction_in
            .filter{it[1][0].type == "protein" || it[1][1].type == "protein"}
            .map{["id": it[0].id]}
        )
    )
    .set{ch_ipsae_in}
    
    IPSAE(
        ch_ipsae_in.map{it[1]},
        ch_ipsae_in.map{it[2]}
    )

    ch_versions = ch_versions.mix(IPSAE.out.versions)

    COLLECT_CONFIDENCE(
        ch_boltz_confidence_scores.ids
        .map{[["id": "all-boltz", "model": "boltz"], it]}
        .mix(ch_colabfold_confidence_scores.ids
            .map{[["id": "all-colabfold", "model": "colabfold"], it]}
        )
        .mix(ch_alphafold3_confidence_scores.ids
            .map{[["id": "all-alphafold3", "model": "alphafold3"], it]}
        ),
        ch_boltz_confidence_scores.json
        .mix(ch_colabfold_confidence_scores.json)
        .mix(ch_alphafold3_confidence_scores.json)
    )
    ch_versions = ch_versions.mix(COLLECT_CONFIDENCE.out.versions)

    
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'proteinfold_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_methods_description   = Channel.value(methodsDescriptionText(ch_multiqc_methods_description))

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files
                        .mix(COLLECT_CONFIDENCE.out.confidence.map{it[1]})
    //ch_multiqc_files.view()
    MULTIQC (
        ch_multiqc_files.collect(sort: true),
        ch_multiqc_config.collect()
            .ifEmpty([]),
        ch_multiqc_custom_config
            .collect()
            .ifEmpty([]),
        ch_multiqc_logo
            .collect()
            .ifEmpty([]),
        [],
        []
    )

    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    emit:
    versions   = ch_versions
    
} 
