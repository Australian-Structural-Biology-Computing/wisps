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

//
// LOCAL MODULE: Boltz
//
include { RUN_BOLTZ } from '../modules/local/run_boltz'
include { RUN_ALPHAFOLD3 } from '../modules/local/run_alphafold3'
include { COLABFOLD_BATCH } from '../modules/local/colabfold_batch'
include { BOLTZ_FASTA } from '../modules/local/data_convertor/boltz_fasta'
include { SPLIT_MSA } from '../modules/local/msa_manager/split_msa'
include { IPSAE } from '../modules/local/ipsae'
include {CREATE_INTERACTIONS}  from '../modules/local/data_convertor/create_interactions'
include { MMSEQS_COLABFOLDSEARCH } from '../modules/local/mmseqs_colabfoldsearch'

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
    ch_samplesheet
    mode
    ch_versions
    ch_boltz_ccd
    ch_boltz_model
    ch_boltz2_aff
    ch_boltz2_conf
    ch_mols
    ch_colabfold_db
    ch_uniref30
    mmseqs_batch_size
    colabfold_model_preset
    ch_colabfold_params
    num_recycles
    ch_alphafold3_params
    tools
    ch_multiqc_config
    ch_multiqc_custom_config
    ch_multiqc_logo
    ch_multiqc_methods_description
    outdir
    analysis_batch_size
    interaction_threshold

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

    ch_samplesheet
    .multiMap{
        if (it[1].extension == "fasta" || it[1].extension == "fa")
            {seq = getFastaSequences(it[1].text).sequence.join(":")}
        else if (it[1].extension == "yaml" || it[1].extension == ".yml")
            {seq = getYamlSequences(it[1].text).sequence.join(":")}
        else if (it[1].extension == "json")
            {seq = ""}

        ids:it[0].id
        types: it[0].type
        groups: it[0].group
        seqs: seq
    }.set{ch_raw_sample_sheet}


    CREATE_INTERACTIONS(
        ch_raw_sample_sheet.ids.collect(),
        ch_raw_sample_sheet.types.collect(),
        ch_raw_sample_sheet.groups.collect(),
        ch_raw_sample_sheet.seqs.collect(),
        interaction_mode,
        interaction_threshold
    )

    ch_samplesheet
    .collect(flat: false)
    .flatMap { data_ls ->
        def list = data_ls.toList()
        def out = []
        for (int i = 0; i < list.size(); i++) {
            for (int j = i; j < list.size(); j++) {
                def key = [list[i][0]['group'].toString(), list[j][0]['group'].toString()].sort().join('-')
                if ((interaction_mode[0] != "all-all" && interaction_mode.contains(key)) || (interaction_mode[0] == "all-all" && (interaction_threshold == 0 || Math.abs(i-j) <= interaction_threshold)))
                {
                    out << [
                                [
                                "id":    list[i][0]["id"]    + "-" + list[j][0]["id"],
                                "group": key
                                ],
                                [list[i][0], list[j][0]],
                                [list[i][1], list[j][1]],
                        ]
                }
            }
        }
        out
    }
    .set{ch_interaction_in}


    ch_interaction_in.count().subscribe{print("total final: ${it}")}

    ch_interaction_in
    .filter{it[1][0].type == "protein" && it[1][1].type == "protein"}
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

    MMSEQS_COLABFOLDSEARCH (
        CREATE_INTERACTIONS.out.interactions.map{[["id": "all_run"], it]},
        ch_colabfold_db,
        ch_uniref30
    )
    ch_versions = ch_versions.mix(MMSEQS_COLABFOLDSEARCH.out.versions)

    MMSEQS_COLABFOLDSEARCH.out.a3m
    .map{it[1]}
    .flatten()
    .map{[it.baseName, it]}
    .cross(
        ch_protein_pairs
        .mix(ch_single_protein)
        .map{[it[0].id, it[0]]}
    )
    .map{[it[1][1], it[0][1]]}
    .set{ch_a3m}
    ch_a3m.first().view()
    ch_protein_pairs.first().view()
    // Prepare interactions for boltz
    ch_boltz_data = Channel.empty()
    ch_split_msa_in = Channel.empty()
    ch_boltz_interactions_in = Channel.empty()

    if ("boltz" in tools.split(",")){
        ch_a3m
        .join(ch_protein_pairs.map{it[0]})
        .set{ch_split_msa_in}
        ch_boltz_interactions_in = ch_interaction_in
    }

    split_batch_id = 0
    SPLIT_MSA(
        ch_split_msa_in
        .map{it[1]}
        .buffer( size: mmseqs_batch_size, remainder: true )
        .map{ split_batch_id += 1; [["id" : "batch-${split_batch_id}-${it.size()}"], it]}
    )

    ch_versions = ch_versions.mix(SPLIT_MSA.out.versions)

    // Adding non protein for boltz
    ch_boltz_data =
        ch_boltz_data
        .mix(ch_boltz_interactions_in
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
        .combine(ch_a3m.filter{it[0].type == "protein"})
        .filter{it[0] == it[2]}
        .map{
            if (it[1][1][0].type == "protein")
                [it[1], [it[3], ""]]
            else
                [it[1], ["", it[3]]]
        }
    )

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
                    it[1] = [it[1][0], it[1][0]]
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

    //prepare interactions for colabfold
    ch_colabfold_interaction_in = Channel.empty()
    if ("colabfold" in tools.split(",")){
        colabfold_batch = 0
        ch_a3m.join(ch_protein_pairs.map{it[0]})
        .map{it[1]}
        .buffer( size: analysis_batch_size, remainder: true )
        .map{
            colabfold_batch += 1;
            [["id": "batch-${colabfold_batch}-${it.size()}"], it]
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


    af3_batch = 0
    ch_alphafold3_interaction_in = Channel.empty()
    if ("alphafold3" in tools.split(",")){
        ch_a3m.join(ch_protein_pairs.map{it[0]})
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


    ipsae_batch = 0

    ch_colabfold_scores
    .join(ch_colabfold_pdb)
    .map{[["id": it[0].id], it]}
    .join(ch_interaction_in
        .filter{it[1][0].type == "protein" || it[1][1].type == "protein"}
        .map{["id": it[0].id]}
    )
    .map{[it[1][1], it[1][2], []]}
    .buffer( size: analysis_batch_size, remainder: true )
    .map{
        ipsae_batch += 1;
        [["id": "batch-${ipsae_batch}-${it.size()}", "model": "colabfold"], it]
    }
    .mix(
        ch_boltz_pae
        .join(ch_boltz_cif)
        .join(ch_boltz_confidence)
        .map{[["id": it[0].id], [it[0], it[1], it[2]], it[3]]}
        .join(ch_interaction_in
            .filter{it[1][0].type == "protein" || it[1][1].type == "protein"}
            .map{["id": it[0].id]}
        )
        .map{[it[1][1], it[1][2], it[2]]}
        .buffer( size: analysis_batch_size, remainder: true )
        .map{
            ipsae_batch += 1;
            [["id": "batch-${ipsae_batch}-${it.size()}", "model": "boltz"], it]
        }
    )
    .mix(
        ch_alphafold3_confidence.join(ch_alphafold3_cif)
        .map{[["id": it[0].id], it]}
        .join(ch_interaction_in
            .filter{it[1][0].type == "protein" || it[1][1].type == "protein"}
            .map{["id": it[0].id]}
        )
        .map{[it[1][1], it[1][2], []]}
        .buffer( size: analysis_batch_size, remainder: true )
        .map{
            ipsae_batch += 1;
            [["id": "batch-${ipsae_batch}-${it.size()}", "model": "alphafold3"], it]
        }
    )
    .set{ch_ipsae_in}

    IPSAE(
        ch_ipsae_in
        .map{[it[0], it[1].collect{[it[0].name, it[1].name]}]},
        ch_ipsae_in
        .map{it[1].flatten()}
    )

    ch_versions = ch_versions.mix(IPSAE.out.versions)

    IPSAE.out.txt
    .transpose()
    .branch {
        boltz: it[0].model == "boltz"
            lines = it[1].text.split("\n").findAll{line -> line.split().size() > 4 && line.split()[4].trim() == "max" }
            max_vals = lines.collect{new BigDecimal(it.split()[5].trim())}
            return [it[1].baseName.split("_")[0], max_vals ? (max_vals.sum() / max_vals.size()) : null]

        alphafold3: it[0].model == "alphafold3"
            lines = it[1].text.split("\n").findAll{line -> line.split().size() > 4 && line.split()[4].trim() == "max" }
            max_vals = lines.collect{new BigDecimal(it.split()[5].trim())}
            return [it[1].baseName.split("_")[0], max_vals ? (max_vals.sum() / max_vals.size()) : null]

        colabfold: it[0].model == "colabfold"
            lines = it[1].text.split("\n").findAll{line -> line.split().size() > 4 && line.split()[4].trim() == "max" }
            max_vals = lines.collect{new BigDecimal(it.split()[5].trim())}
            return [it[1].baseName.split("_")[0], max_vals ? (max_vals.sum() / max_vals.size()) : null]
    }
    .set{ch_ipsae_out}

    ch_interaction_in.map{it[0].id}
    .join(ch_ipsae_out.boltz, remainder: true)
    .join(ch_ipsae_out.colabfold, remainder: true)
    .join(ch_ipsae_out.alphafold3, remainder: true)
    .map{"${it[0]},${it[1] != null ? it[1] : ""},${ it[2] != null ? it[2] : ""},${ it[3] != null ? it[3] : ""}\n"}
        .toSortedList()
        .flatten()
        .collectFile(
            storeDir: "${outdir}/ipsae",
            name: 'ipsae_scores.csv',
            seed: "Sample,boltz_ipsae,colabfold_ipsae,alphafold3_ipsae\n")
    .set{ch_ipsae_scores}

    ch_boltz_confidence
    .collect(flat: false, sort: true).multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_boltz_confidence_scores}

    ch_alphafold3_summary_confidences
    .collect(flat: false, sort: true).multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_alphafold3_confidence_scores}

    ch_colabfold_scores
    .collect(flat: false, sort: true).multiMap{json_list ->
        ids: json_list.collect{it[0].id}
        json: json_list.collect{it[1]}
    }.set{ch_colabfold_confidence_scores}


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
                        .mix(ch_ipsae_scores)
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
