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
include { MSA } from '../subworkflows/local/msa'

//
// MODULE: Boltz
//
include { RUN_BOLTZ } from '../modules/local/run_boltz'

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
    ch_colabfold_db // channel: [ path(colabfold_db) ]
    ch_uniref30     // channel: [ path(uniref30) ]
    mmseq_batch_size // number

    main:
    ch_multiqc_files = Channel.empty()
    
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
    .branch {
        fasta: it[1].extension == "fasta" || it[1].extension == "fa"
            it[0].cnt = getFastaSequences(it[1].text).size()
            return it
        yaml: it[1].extension == "yaml" || it[1].extension == ".yml"
            it[0].cnt = getYamlSequences(it[1].text).size()
            return it
        json: it[1].extension == "json"
    }
    .set{ch_input}
    

    MSA (
        ch_input.fasta
        .mix(ch_input.json)
        .mix(ch_input.yaml)
        .filter{it[0].type == "protein"},
        ch_colabfold_db,
        ch_uniref30,
        mmseq_batch_size
    )
    ch_versions = ch_versions.mix(MSA.out.versions)
    
    ch_input.fasta
        .mix(ch_input.json)
        .mix(ch_input.yaml)
        .map{[it[0].id, it]}
        .join(
            MSA.out.a3m.map{[it[0].id, it]}
            , remainder: true
        )
        .map{[it[2][0], it[1][1], it[2][1]]}
        .set{ch_input_aligned}
        
    ch_input_aligned
    .combine(ch_input_aligned)
    .map{[
            [  
                "id": [it[0]["id"], it[3]["id"]].min() + "-" + [it[0]["id"], it[3]["id"]].max(),
                "group": [it[0]["group"], it[3]["group"]].min() + "-" + [it[0]["group"], it[3]["group"]].max(),
            ],
            [it[1].name, it[4].name], 
            [it[2] ? it[2].name : "", it[5] ? it[5].name : ""], 
            [it[0]["type"], it[3]["type"]],
            [it[1], it[2], it[4], it[5]].findAll{it}.unique{it.toUriString()},

    ]}.set{ch_pairs}
    
    //ch_pairs.view()
    
    ch_pairs
    .join(ch_unique_pairs)
    .set{ch_interaction_input}
    
    //ch_interaction_input.view()

    PREPARE_INTERACTIONS(
        ch_interaction_input.map{[it[0], it[1]]},
        ch_interaction_input.map{[it[0], it[2]]},
        ch_interaction_input.map{[it[0], it[3]]},
        ch_interaction_input.map{it[4]}
    )
    
    PREPARE_INTERACTIONS.out.fasta
    .join(ch_interaction_input)
    .set{ch_boltz_input}
    
    RUN_BOLTZ(
        ch_boltz_input.map{[it[0], it[1]]},
        ch_boltz_input.map{it[5].findAll { it.name.endsWith(".a3m") }},
        ch_boltz_model,
        ch_boltz_ccd
    )
    RUN_BOLTZ.out.confidence.toSortedList().multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_confidence_scores}

    COLLECT_CONFIDENCE(
        ch_confidence_scores.ids.map{[["id": "full_run"], it]},
        ch_confidence_scores.json
    )
    
    //ch_confidence_scores.ids.view()
    //ch_confidence_scores.json.view()

    emit:
    versions   = ch_versions
    msa        = RUN_BOLTZ.out.msa
    structures = RUN_BOLTZ.out.structures
    confidence = RUN_BOLTZ.out.confidence
    run_confidence = COLLECT_CONFIDENCE.out.confidence
    plddt      = RUN_BOLTZ.out.plddt
    pdb        = RUN_BOLTZ.out.pdb

} 
