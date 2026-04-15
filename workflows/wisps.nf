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
include { COLLECT_CONFIDENCE as COLLECT_CONFIDENCE_BOLTZ } from '../modules/local/collect_confidence'
include { COLLECT_CONFIDENCE as COLLECT_CONFIDENCE_COLABFOLD } from '../modules/local/collect_confidence'
include { COLLECT_CONFIDENCE as COLLECT_CONFIDENCE_AF3 } from '../modules/local/collect_confidence'
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
include { IPSAE } from '../modules/local/ipsae'
include {CREATE_INTERACTIONS}  from '../modules/local/data_convertor/create_interactions'
include {CREATE_INTERACTION_POOLS} from '../modules/local/data_convertor/create_interaction_pools'
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
    ch_boltz2_aff
    ch_boltz2_conf
    ch_mols
    ch_colabfold_envdb
    ch_colabfold_uniref30
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
    colabfold_batch_size
    interaction_neighbours
    pool
    pool_size

    main:
    ch_multiqc_files = Channel.empty()
    ch_confidence_scores = Channel.empty()

    use_interaction_pools = (pool instanceof Boolean) ? pool : pool.toString().toBoolean()
    interaction_mode = mode.split(",").collect { pair ->
        pair.split('-').sort().join('-')
    }
    if ("all-all" in interaction_mode) {
        interaction_mode = ["all-all"]
    }

    if (use_interaction_pools) {
        if ((pool_size as Integer) <= 0) {
            error("When using --pool, --pool_size must be a positive integer. Received: ${pool_size}")
        }
        if ("manual" in interaction_mode || "all-all" in interaction_mode) {
            error("When using --pool, --mode must be explicit group pairs (e.g. 'A-B' or 'A-A,A-B'). 'manual' and 'all-all' are not supported.")
        }

        required_pool_groups = interaction_mode.collect{it.split("-")}.flatten().unique()
        ch_samplesheet
            .map{it[0].group}
            .unique()
            .toList()
            .subscribe{
                emitted -> {
                    if (required_pool_groups - emitted){
                        log.error "The groups ${required_pool_groups - emitted} are not in the sample sheet!"
                        exit 1
                    }
                }
            }
        log.info "Using interaction pooling with modes: ${interaction_mode.join(',')} and max total sequence length ${pool_size}"
    } else {
        if (!("manual" in interaction_mode) && !("all-all" in interaction_mode)) {
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
    }

    ch_samplesheet
    .multiMap{
        if (it[1].extension == "fasta" || it[1].extension == "fa")
            {seq = getFastaSequences(it[1].text).sequence.join(":")}
        else if (it[1].extension == "yaml" || it[1].extension == "yml")
            {seq = getYamlSequences(it[1].text).sequence.join(":")}
        else if (it[1].extension == "json")
            {seq = ""}

        ids:it[0].id
        types: it[0].type
        groups: it[0].group
        seqs: seq
    }.set{ch_raw_sample_sheet}


    ch_interaction_raw = Channel.empty()
    ch_interaction_info = Channel.empty()

    if (use_interaction_pools) {
        CREATE_INTERACTION_POOLS(
            ch_raw_sample_sheet.ids.collect(),
            ch_raw_sample_sheet.types.collect(),
            ch_raw_sample_sheet.groups.collect(),
            ch_raw_sample_sheet.seqs.collect(),
            interaction_mode,
            pool_size
        )

        ch_interaction_info = CREATE_INTERACTION_POOLS.out.interaction_mapping
        .flatMap { mapping_tsv ->
            def out = []
            def lines = mapping_tsv.readLines().findAll { it?.trim() }
            if (lines.isEmpty()) {
                throw new IllegalStateException("Invalid interaction mapping file ${mapping_tsv}: file is empty.")
            }

            lines.drop(1).each { line ->
                def cols = line.split("\t", -1)
                if (cols.size() < 3) {
                    throw new IllegalStateException("Invalid interaction mapping row in ${mapping_tsv}: '${line}'")
                }
                def interaction_id = cols[0].trim()
                def int_chain_map = cols[1].trim()
                def str_chain_map = cols[2].trim()
                def report_id = cols.size() > 3 && cols[3].trim() ? cols[3].trim() : interaction_id
                def left_source_id = cols.size() > 4 ? cols[4].trim() : ""
                def right_source_id = cols.size() > 5 ? cols[5].trim() : ""
                out << [interaction_id, int_chain_map, str_chain_map, report_id, left_source_id, right_source_id]
            }
            out
        }

        ch_interaction_raw = CREATE_INTERACTION_POOLS.out.interactions
    } else if (!("manual" in interaction_mode)) {
        CREATE_INTERACTIONS(
            ch_raw_sample_sheet.ids.collect(),
            ch_raw_sample_sheet.types.collect(),
            ch_raw_sample_sheet.groups.collect(),
            ch_raw_sample_sheet.seqs.collect(),
            interaction_mode,
            interaction_neighbours
        )

        ch_interaction_info = CREATE_INTERACTIONS.out.interaction_mapping
        .flatMap { mapping_tsv ->
            def out = []
            def lines = mapping_tsv.readLines().findAll { it?.trim() }
            if (lines.isEmpty()) {
                throw new IllegalStateException("Invalid interaction mapping file ${mapping_tsv}: file is empty.")
            }

            lines.drop(1).each { line ->
                def cols = line.split("\t", -1)
                if (cols.size() != 3) {
                    throw new IllegalStateException("Invalid interaction mapping row in ${mapping_tsv}: '${line}'")
                }
                def interaction_id = cols[0].trim()
                def int_chain_map = cols[1].trim()
                def str_chain_map = cols[2].trim()
                out << [interaction_id, int_chain_map, str_chain_map, interaction_id, "", ""]
            }
            out
        }

        ch_interaction_raw = CREATE_INTERACTIONS.out.interactions
    } else {
        ch_interaction_raw = ch_samplesheet
            .map { it[1] }
            .filter { it.extension == "fasta" || it.extension == "fa" }
            .ifEmpty { error("Manual mode requires at least one FASTA file in the samplesheet sequence column.") }
            .collectFile(name: "manual_interactions.fasta", newLine: false) { fasta_file ->
                ["manual_interactions.fasta", fasta_file.text.trim() + "\n"]
            }
    }

    ch_interaction_has_protein = ch_interaction_raw.flatMap { interactions_fasta ->
        def out = []
        def lines = interactions_fasta.readLines().findAll { it?.trim() }
        if (lines.size() % 2 != 0) {
            throw new IllegalStateException("Invalid interactions FASTA layout in ${interactions_fasta}: expected an even number of non-empty lines.")
        }

        lines.collate(2).each { chunk ->
            def header = chunk[0]
            if (!header.startsWith(">")) {
                throw new IllegalStateException("Invalid FASTA header in ${interactions_fasta}: '${header}'")
            }

            def interaction_id = header.substring(1).trim()
            def sequence_line = chunk[1].trim()
            def has_protein = sequence_line.split(":").count { !it.contains("|") } > 1
            out << [interaction_id, has_protein]
        }
        out
    }

    if ("manual" in interaction_mode) {
        ch_interaction_info = ch_interaction_has_protein
            .map { interaction_id, _ -> [interaction_id, "", "", interaction_id, "", ""] }
    }

    ch_interaction_has_protein_map = ch_interaction_has_protein
        .collect(flat: false)
        .map { rows ->
            rows.collectEntries { row -> [(row[0]): row[1]] }
        }

    ch_interaction_in = ch_interaction_info
        .combine(ch_interaction_has_protein_map)
        .map { interaction_id, int_chain_map, str_chain_map, report_id, left_source_id, right_source_id, has_protein_map ->
            def has_protein = has_protein_map.containsKey(interaction_id) ? has_protein_map[interaction_id] : false
            [[
                "id": interaction_id,
                "report_id": report_id,
                "int_chain_map": int_chain_map,
                "str_chain_map": str_chain_map,
                "left_source_id": left_source_id,
                "right_source_id": right_source_id
            ], has_protein]
        }

    ch_interaction_in.count().subscribe{print("total final: ${it}")}

    MMSEQS_COLABFOLDSEARCH (
        ch_interaction_raw.map{[["id": "all_run"], it]},
        ch_colabfold_envdb,
        ch_colabfold_uniref30
    )
    ch_versions = ch_versions.mix(MMSEQS_COLABFOLDSEARCH.out.versions)

    MMSEQS_COLABFOLDSEARCH.out.a3m
    .map{it[1]}
    .flatten()
    .map{[it.baseName, it]}
    .set{ch_a3m}

    MMSEQS_COLABFOLDSEARCH.out.json
    .map{it[1]}
    .flatten()
    .map{[it.baseName, it]}
    .set{ch_af3_json}


    // Prepare interactions for boltz directly from MMSEQS outputs
    ch_boltz_in = Channel.empty()
    if ("boltz" in tools.split(",")){
        MMSEQS_COLABFOLDSEARCH.out.yaml
        .map{it[1]}
        .flatten()
        .map{[it.baseName, it]}
        .join(
            MMSEQS_COLABFOLDSEARCH.out.msa_csv
            .map{it[1]}
            .flatten()
            .map{[it.baseName.replaceFirst(/_[A-Z]+$/, '').replace("_boltz_interaction_input", ""), it]}
            .groupTuple(sort: true),
            remainder: true
        )
        .map{id, yaml, msa_csv -> [id, yaml, msa_csv ?: []]}
        .buffer( size: analysis_batch_size, remainder: true )
        .set{ch_boltz_in}
    }

    boltz_batch = 0

    RUN_BOLTZ(
        ch_boltz_in.map{boltz_batch += 1; [["id": "batch-${boltz_batch}-${it.size()}"], it.collect{it[1]}]},
        ch_boltz_in.map{it.collect{it[2]}.flatten().findAll{it}.unique{it.toUriString()}},
        ch_boltz2_aff,
        ch_boltz2_conf,
        ch_mols
    )


    RUN_BOLTZ.out.confidence
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_model_0")[0].split("confidence_")[1], "model": "boltz"], it]}
    .set{ch_boltz_confidence}
    
    RUN_BOLTZ.out.affinity
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_model_0")[0].split("affinity_")[1], "model": "boltz"], it]}
    .set{ch_boltz_affinity}

    RUN_BOLTZ.out.pae
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_model_0")[0].split("pae_")[1], "model": "boltz"], it]}
    .set{ch_boltz_pae}

    RUN_BOLTZ.out.cif
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_model_0")[0], "model": "boltz"], it]}
    .set{ch_boltz_cif}

    //prepare interactions for colabfold
    ch_colabfold_interaction_in = Channel.empty()
    if ("colabfold" in tools.split(",")){
        colabfold_batch = 0
        //ch_a3m.join(ch_protein_pairs.map{it[0]})
        ch_a3m
        .map{it[1]}
        .buffer( size: colabfold_batch_size, remainder: true )
        .map{
            colabfold_batch += 1;
            [["id": "batch-${colabfold_batch}-${it.size()}"], it]
        }
        .set{
            ch_colabfold_interaction_in
        }
    }

    COLABFOLD_BATCH(
        ch_colabfold_interaction_in,
        ch_colabfold_params,
        num_recycles
    )
    ch_versions = ch_versions.mix(COLABFOLD_BATCH.out.versions)

    COLABFOLD_BATCH.out.top_ranked_scores
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_scores_rank_")[0], "model": "colabfold"], it]}
    .set{ch_colabfold_scores}

    COLABFOLD_BATCH.out.pdb
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_unrelaxed_")[0], "model": "colabfold"], it]}
    .set{ch_colabfold_pdb}


    af3_batch = 0
    ch_alphafold3_interaction_in = Channel.empty()
    if ("alphafold3" in tools.split(",")){
        //ch_af3_json.join(ch_protein_pairs.map{it[0]})
        ch_af3_json
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
    //since *_confidences overmatches outputs
    ch_jsons = RUN_ALPHAFOLD3.out.jsons
        .flatMap { meta, files ->
            files.collect { f -> tuple(meta, f) }
        }

    confidence = ch_jsons
        .filter { meta, f ->
            f.name.endsWith('_confidences.json') &&
            !f.name.endsWith('_summary_confidences.json')
        }

    summary_confidence = ch_jsons
        .filter { meta, f ->
            f.name.endsWith('_summary_confidences.json')
        }

    RUN_ALPHAFOLD3.out.top_ranked_cif
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_model")[0], "model": "alphafold3"], it]}
    .set{ch_alphafold3_cif}

    confidence
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_confidences")[0], "model": "alphafold3"], it]}
    .set{ch_alphafold3_confidence}

    summary_confidence
    .map{it[1]}
    .flatten()
    .map{[["id": it.baseName.split("_summary_confidences")[0], "model": "alphafold3"], it]}
    .set{ch_alphafold3_summary_confidences}


    ipsae_batch = 0
    ch_protein_interaction_by_id = ch_interaction_has_protein
        .filter { it[1] }
        .map { [it[0], true] }

    ch_colabfold_scores
    .join(ch_colabfold_pdb)
    .map { [it[0].id, it[1], it[2]] }
    .join(ch_protein_interaction_by_id)
    .map { id, score_json, pdb_file, _ -> [score_json, pdb_file, []] }
    .buffer( size: analysis_batch_size, remainder: true )
    .map{
        ipsae_batch += 1;
        [["id": "batch-${ipsae_batch}-${it.size()}", "model": "colabfold"], it]
    }
    .mix(
        ch_boltz_pae
        .join(ch_boltz_cif)
        .join(ch_boltz_confidence)
        .map { [it[0].id, it[1], it[2], it[3]] }
        .join(ch_protein_interaction_by_id)
        .map { id, pae_file, cif_file, confidence_json, _ -> [pae_file, cif_file, confidence_json] }
        .buffer( size: analysis_batch_size, remainder: true )
        .map{
            ipsae_batch += 1;
            [["id": "batch-${ipsae_batch}-${it.size()}", "model": "boltz"], it]
        }
    )
    .mix(
        ch_alphafold3_confidence.join(ch_alphafold3_cif)
        .map { [it[0].id, it[1], it[2]] }
        .join(ch_protein_interaction_by_id)
        .map { id, confidence_json, cif_file, _ -> [confidence_json, cif_file, []] }
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

    def parseIpsaeMaxRows = { file, model, sampleId ->
        def lines = file.text
            .split("\n")
            .findAll { line ->
                def cols = line.split()
                cols.size() > 5 && cols[4].trim() == "max"
            }
        def entries = lines.collect { line ->
            def cols = line.split()
            def pairId = "${cols[0].trim()}_${cols[1].trim()}"
            [pairId, new BigDecimal(cols[5].trim())]
        }
        [sampleId, entries]
    }

    def parseStrChainMap = { strChainMap ->
        if (!strChainMap || !strChainMap.contains(":")) {
            return [[], []]
        }
        def parts = strChainMap.split(":", 2)
        def leftIds = parts[0].split(",").collect { it.trim() }.findAll { it != "" }
        def rightIds = parts[1].split(",").collect { it.trim() }.findAll { it != "" }
        [leftIds, rightIds]
    }

    def computeMaxCrossGroupIpsae = { entries, strChainMap ->
        if (!entries) {
            return null
        }
        def overallMax = entries.collect { it[1] }.findAll { it != null }.max()

        // If no chain map exists, fall back to the overall max IPSAE.
        if (!strChainMap || !strChainMap.contains(":")) {
            return overallMax
        }

        def parsed = parseStrChainMap(strChainMap)
        def leftIds = parsed[0] as Set
        def rightIds = parsed[1] as Set
        if (leftIds.isEmpty() || rightIds.isEmpty()) {
            return overallMax
        }

        def crossScores = entries.collect { pair ->
            def pairId = pair[0].toString()
            def score = pair[1]
            def chains = pairId.split("_", 2)
            if (chains.size() != 2) {
                chains = pairId.split("-", 2)
            }
            if (chains.size() != 2) {
                return null
            }
            def left = chains[0].trim()
            def right = chains[1].trim()
            ((leftIds.contains(left) && rightIds.contains(right)) ||
            (leftIds.contains(right) && rightIds.contains(left))) ? score : null
        }.findAll { it != null }

        crossScores ? crossScores.max() : null
    }

    IPSAE.out.txt
    .transpose()
    .branch {
        boltz: it[0].model == "boltz"
            return parseIpsaeMaxRows(it[1], "boltz", it[1].baseName.split("_model_0")[0])

        alphafold3: it[0].model == "alphafold3"
            return parseIpsaeMaxRows(it[1], "alphafold3", it[1].baseName.split("_model_")[0])

        colabfold: it[0].model == "colabfold"
            return parseIpsaeMaxRows(it[1], "colabfold", it[1].baseName.split("_unrelaxed_")[0])
    }
    .set{ch_ipsae_out}

    def active_modes = tools
    .split(',')
    .collect { it.trim() }
    .findAll { it in ['boltz','colabfold','alphafold3'] }

    // 1) Long format: [sample_id, pair_id, model, score]
    ch_report_meta_grouped_by_input_id = ch_interaction_in
        .map { [it[0].id, it[0]] }
        .groupTuple()
        .map { inputId, metas -> [inputId, metas] }

    ch_long = Channel.empty()
    if ('boltz' in active_modes)
        ch_long = ch_long.mix(ch_ipsae_out.boltz.flatMap { id, entries -> entries.collect { e -> [id, e[0], 'boltz', e[1]] } })
    if ('colabfold' in active_modes)
        ch_long = ch_long.mix(ch_ipsae_out.colabfold.flatMap { id, entries -> entries.collect { e -> [id, e[0], 'colabfold', e[1]] } })
    if ('alphafold3' in active_modes)
        ch_long = ch_long.mix(ch_ipsae_out.alphafold3.flatMap { id, entries -> entries.collect { e -> [id, e[0], 'alphafold3', e[1]] } })
    ch_long = ch_long
        .join(ch_report_meta_grouped_by_input_id, remainder: true)
        .flatMap { inputId, pairId, model, score, metaList ->
            def metas = (metaList instanceof List && !metaList.isEmpty())
                ? metaList
                : [[report_id: inputId, str_chain_map: ""]]
            metas.collect { interaction_meta ->
                def reportId = interaction_meta?.report_id ?: inputId
                [reportId, pairId, model, score]
            }
        }

    // 2) Pivot per pair: [pair_id, sampleModelScoreMap]
    ch_pair_wide = ch_long
        .map { sampleId, pairId, model, score -> [pairId, [sampleId, model, score]] }
        .groupTuple()
        .map { pairId, rows ->
            def bySample = [:].withDefault { [:] }
            rows.each { r ->
                bySample[r[0]][r[1]] = r[2]
            }
            [pairId, bySample]
        }

    // 3) Create one CSV per pair with model-only headers
    ch_interaction_in.map { it[0].report_id }.unique().toSortedList().map { [sample_ids: it] }
        .combine(ch_pair_wide.collect(flat: false).map { [pairs: it] })
        .flatMap { left, right ->
            def ids = left.sample_ids
            def pairs = right.pairs
            pairs.collect { pairEntry ->
                def pairId = pairEntry[0]
                def bySample = pairEntry[1]
                def header = "Sample," + active_modes.join(",") + "\n"
                def body = ids.collect { id ->
                    def mm = bySample[id] ?: [:]
                    def vals = active_modes.collect { mode -> mm[mode] != null ? mm[mode] : "" }
                    "${id},${vals.join(',')}\n"
                }.join("")
                def safePairId = pairId.replaceAll(/[^A-Za-z0-9._-]+/, "_")
                ["ipsae_scores_${safePairId}.csv", header + body]
            }
        }
        .collectFile(
            storeDir: "${outdir}/ipsae",
            newLine: false
        ) { item -> [item[0], item[1]] }
        .set { ch_ipsae_scores_by_pair }

    // 4) Additional summary file: max cross-group (left/right in str_chain_map) IPSAE per sample/model
    ch_ipsae_entries = Channel.empty()
    if ('boltz' in active_modes)
        ch_ipsae_entries = ch_ipsae_entries.mix(ch_ipsae_out.boltz.map { inputId, entries -> [inputId, "boltz", entries] })
    if ('colabfold' in active_modes)
        ch_ipsae_entries = ch_ipsae_entries.mix(ch_ipsae_out.colabfold.map { inputId, entries -> [inputId, "colabfold", entries] })
    if ('alphafold3' in active_modes)
        ch_ipsae_entries = ch_ipsae_entries.mix(ch_ipsae_out.alphafold3.map { inputId, entries -> [inputId, "alphafold3", entries] })

    ch_model_rows = ch_ipsae_entries
        .combine(ch_report_meta_grouped_by_input_id, by: 0)
        .flatMap { inputId, model, entries, metaList ->
            def metas = (metaList instanceof List && !metaList.isEmpty())
                ? metaList
                : [[report_id: inputId, str_chain_map: ""]]
            metas.collect { interaction_meta ->
                def reportId = interaction_meta?.report_id ?: inputId
                def strChainMap = interaction_meta?.str_chain_map ?: ""
                [[reportId, model], computeMaxCrossGroupIpsae(entries, strChainMap)]
            }
        }

    ch_model_wide = ch_model_rows
        .groupTuple()
        .map { key, scores -> [key[0], [key[1], scores.max()]] }
        .groupTuple()
        .map { sampleId, pairs ->
            def m = pairs.collectEntries { p -> [(p[0]): p[1]] }
            [sampleId, m]
        }

    def summary_header = "Sample," + active_modes.join(",") + "\n"
    ch_interaction_in.map { it[0].report_id }.unique().map { [it, true] }
        .join(ch_model_wide, remainder: true)
        .map { id, _, m ->
            def mm = m ?: [:]
            def vals = active_modes.collect { mode -> mm[mode] != null ? mm[mode] : "" }
            "${id},${vals.join(',')}\n"
        }
        .toSortedList()
        .flatMap { it }
        .collectFile(
            storeDir: "${outdir}/ipsae",
            name: 'ipsae_scores.csv',
            seed: summary_header
        )
        .set { ch_ipsae_scores_max }

    ch_ipsae_scores = ch_ipsae_scores_by_pair.mix(ch_ipsae_scores_max)

    ch_confidence_meta = ch_interaction_in
        .map { it[0] }
        .unique { it.report_id ?: "${it.id}:${it.int_chain_map}:${it.str_chain_map}" }
        .collect(flat: false, sort: true)

    ch_boltz_confidence_json = ch_boltz_confidence
        .map { it[1] }
        .mix(ch_boltz_affinity.map { it[1] })
        .collect(flat: false, sort: true)

    ch_alphafold3_confidence_json = ch_alphafold3_summary_confidences
        .map { it[1] }
        .collect(flat: false, sort: true)

    ch_colabfold_confidence_json = ch_colabfold_scores
        .map { it[1] }
        .collect(flat: false, sort: true)


    COLLECT_CONFIDENCE_BOLTZ(
        ch_confidence_meta.map{[["id": "all-boltz", "model": "boltz"], it]},
        ch_boltz_confidence_json
    )
    ch_versions = ch_versions.mix(COLLECT_CONFIDENCE_BOLTZ.out.versions)
    ch_confidence_scores_all = COLLECT_CONFIDENCE_BOLTZ.out.confidence

    COLLECT_CONFIDENCE_COLABFOLD(
        ch_confidence_meta.map{[["id": "all-colabfold", "model": "colabfold"], it]},
        ch_colabfold_confidence_json
    )
    ch_versions = ch_versions.mix(COLLECT_CONFIDENCE_COLABFOLD.out.versions)
    ch_confidence_scores_all = ch_confidence_scores_all.mix(COLLECT_CONFIDENCE_COLABFOLD.out.confidence)

    COLLECT_CONFIDENCE_AF3(
        ch_confidence_meta.map{[["id": "all-alphafold3", "model": "alphafold3"], it]},
        ch_alphafold3_confidence_json
    )
    ch_versions = ch_versions.mix(COLLECT_CONFIDENCE_AF3.out.versions)
    ch_confidence_scores_all = ch_confidence_scores_all.mix(COLLECT_CONFIDENCE_AF3.out.confidence)



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
                        .mix(ch_confidence_scores_all.map{it[1]})
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
