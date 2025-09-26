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
    ch_colabfold_db // channel: [ path(colabfold_db) ]
    ch_boltz2_aff   // channel: [ path(boltz2_aff) ]
    ch_boltz2_conf  // channel: [ path(boltz2_conf) ]
    ch_mols         // channel: [ path(mols) ]
    ch_uniref30     // channel: [ path(uniref30) ]
    mmseqs_batch_size // number
    colabfold_model_preset
    ch_colabfold_params
    num_recycles
    ch_alphafold3_params // channel: path(alphafold2_params)
    ch_small_bfd         // channel: path(small_bfd)
    ch_mgnify            // channel: path(mgnify)
    ch_mmcif_files       // channel: path(mmcif_files)
    ch_uniref90          // channel: path(uniref90)
    ch_pdb_seqres        // channel: path(pdb_seqres)
    ch_uniprot           // channel: path(uniprot)
    tools

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

    /*ch_samplesheet
    .map {
        if (it[1].extension == "fasta" || it[1].extension == "fa")
            {it[0].cnt = getFastaSequences(it[1].text).size()}
        else if (it[1].extension == "yaml" || it[1].extension == ".yml")
            {it[0].cnt = getYamlSequences(it[1].text).size()}
        else if (it[1].extension == "json")
            {it[0].cnt = 0}
    }
    .set{ch_input_seqs}
    */
    
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
    //ch_protein_pairs.view()
    MSA (
        ch_protein_pairs
        .mix(
            ch_interaction_in
            .map{
                if (it[1][0].type != "protein" && it[1][1].type == "protein"){
                    [it[1][1], it[2][1]]   
                }else if (it[1][1].type != "protein" && it[1][0].type == "protein"){
                    [it[1][0], it[2][0]]   
                }
            }.unique()
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

    ch_boltz_data = ch_boltz_data.mix(ch_boltz_interactions_in
        .filter{it[1][0].type != "protein" && it[1][1].type != "protein"}
        .map{[it, ["", ""]]}
    )

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

    SPLIT_MSA(
        ch_split_msa_in.join(ch_protein_pairs.map{it[0]})
    )

    ch_versions = ch_versions.mix(SPLIT_MSA.out.versions)
    
    ch_boltz_data = ch_boltz_data.mix(ch_boltz_interactions_in
        .filter{it[1][0].type == "protein" && it[1][1].type == "protein"}
        .map{[["id": it[0].id], it]}
        .join(
            SPLIT_MSA.out.msa_csv
            .map{
                if (it[1] instanceof List){
                    if (it[1].size() != 2){
                        log.warn "Something wrong with the data expected 2 csv!"
                    }
                }else{
                    it[1] = [it[1], it[1]] 
                }
                it
            }
            
        )
        .map{[it[1], it[2]]}
    )
    
    PREPARE_INTERACTIONS(
        ch_boltz_data.map{[it[0][0], it[0][2].collect{it.name}]},
        ch_boltz_data.map{[it[0][0], it[1].collect{it ? it.name: ""}]},
        ch_boltz_data.map{[it[0][0], [it[0][1][0].type, it[0][1][1].type]]},
        ch_boltz_data.map{[it[0][2][0], it[0][2][1], it[1][0], it[1][1]].findAll{it}.unique{it.toUriString()}}
    )
    //PREPARE_INTERACTIONS.out.fasta.view()

    PREPARE_INTERACTIONS.out.fasta.join(
        ch_boltz_data.map{[it[0][0], [it[1][0], it[1][1]].findAll{it}.unique{it.toUriString()}]}
    )
    .map{
            meta = it[0].clone()
            meta.model = "boltz"
            [meta, it[1], it[2]]
    }
    .set{ch_boltz_in}

    //ch_boltz_in.view()
    
    
    RUN_BOLTZ(
        ch_boltz_in.map{[it[0], it[1]]},
        ch_boltz_in.map{it[2]},
        ch_boltz_model,
        ch_boltz_ccd,
        ch_boltz2_aff,
        ch_boltz2_conf,
        ch_mols
    )
    
    RUN_BOLTZ.out.confidence.collect(flat: false, sort: true).multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_confidence_scores}
    
    //prepare interactions for colabfold
    ch_colabfold_interaction_in = Channel.empty()
    if ("colabfold" in tools.split(",")){
        MSA.out.a3m.join(ch_protein_pairs.map{it[0]})
        .map{
            meta = it[0].clone()
            meta.model = "colabfold"
            [meta, it[1]]
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
    
    ch_alphafold3_interaction_in = Channel.empty()
    if ("ZZalphafold3" in tools.split(",")){
        MSA.out.json.join(ch_protein_pairs.map{it[0]})
        .map{
            meta = it[0].clone()
            meta.model = "alphafold3"
            [meta, it[1]]}
       .set{
            ch_alphafold3_interaction_in
        }
    }


    RUN_ALPHAFOLD3 (
        ch_alphafold3_interaction_in,
        ch_alphafold3_params,
        ch_small_bfd,
        ch_mgnify,
        ch_mmcif_files,
        ch_uniref90,
        ch_pdb_seqres,
        ch_uniprot
    )
    
    IPSAE(
        COLABFOLD_BATCH.out.top_ranked_scores.join(
            COLABFOLD_BATCH.out.top_ranked_pdb
        ).map{[["id": it[0].id], it]}
            .join(ch_protein_pairs.map{it[0]})
        .mix(
            RUN_BOLTZ.out.pae.join(RUN_BOLTZ.out.cif)
            .map{[["id": it[0].id], it]}
            .join(ch_protein_pairs.map{it[0]})
        )
        .mix(
            RUN_ALPHAFOLD3.out.pae.join(RUN_ALPHAFOLD3.out.top_ranked_cif)
            .map{[["id": it[0].id], it]}
            .join(ch_protein_pairs.map{it[0]})
        )
        .map{it[1]}    
    )


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
