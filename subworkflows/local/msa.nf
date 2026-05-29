//
// Post processing analysis for the predicted structures
//

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { MMSEQS_COLABFOLDSEARCH } from '../../modules/local/mmseqs_colabfoldsearch'

workflow MSA {

    take:
    ch_samplesheet
    ch_colabfold_envdb        // channel: path(colabfold_envdb)
    ch_colabfold_uniref30            // channel: path(colabfold_uniref30)
    batch_size

    main:
    ch_versions = Channel.empty()
    ch_a3m      = Channel.empty()
    ch_json      = Channel.empty()

    ch_samplesheet
    .branch {
        fasta: it[1].extension == "fasta" || it[1].extension == "fa"
        yaml: it[1].extension == "yaml" || it[1].extension == ".yml"
        json: it[1].extension == "json"
    }
    .set{ch_input}

    ch_input.yaml
    .map{
        [it[0], "${it[0].id},${getYamlSequences(it[1].text).collect { it.sequence }.join(':')}"]
    }
    .set{ch_yaml_seqs}

    if (batch_size > 1){
        def batch_itr = 0

        ch_input.fasta
        .map{
            "${it[0].id},${getFastaSequences(it[1].text).collect { it.sequence }.join(':')}"
        }
        .mix(
            ch_yaml_seqs.map{it[1]}
        )
        .buffer( size: batch_size, remainder: true )
        .collectFile {
            batch_itr += 1;
            [ "input_seqs_${batch_itr}.csv", "id,sequence\n" + it.join("\n") + '\n' ]
        }
        .map{[["id": it.baseName], it]}
        .set {ch_input_seqs}
    }else{
        ch_input.fasta
        .mix(
            ch_yaml_seqs
            .collectFile {
                [ "${it[0].id}.csv", "id,sequence\n" + it[1].join("\n") + '\n' ]
            }
            .map{[["id": it.baseName], it]}
            .join(
                ch_yaml_seqs.map{[it[0].id, it[0]]}
            )
            .map{[it[2], it[1]]}
        ).set{ch_input_seqs}
    }
    //ch_input_seqs.view()

    MMSEQS_COLABFOLDSEARCH (
        ch_input_seqs,
        ch_colabfold_envdb,
        ch_colabfold_uniref30
    )
    ch_versions = ch_versions.mix(MMSEQS_COLABFOLDSEARCH.out.versions)

    ch_a3m = ch_a3m.mix(
        ch_input.fasta.mix(ch_input.yaml)
        .map{[it[1].baseName, it[0]]}
        .combine(
            MMSEQS_COLABFOLDSEARCH.out.a3m
            .map{it[1]}
            .flatten()
            .map {[it.baseName, it]},
            by:0
        )
        .map{[it[1], it[2]]}
    )
    ch_json = ch_json.mix(
        ch_input.fasta.mix(ch_input.yaml)
        .map{[it[1].baseName, it[0]]}
        .combine(
            MMSEQS_COLABFOLDSEARCH.out.json
            .map{it[1]}
            .flatten()
            .map {[it.baseName, it]},
            by:0
        )
        .map{[it[1], it[2]]}
    )

    emit:
    formated_input          = ch_input.fasta.mix(ch_input.yaml)
    a3m            = ch_a3m
    json    = ch_json
    versions       = ch_versions
}

def getYamlSequences(yamlData) {
    List<Map> enrichedEntries = []
    Map currentEntry = [:]
    inSequences = false
    yamlData.split("\n").each { line ->
        def trimmed = line.trim()

        // Detect start of sequences section
        if (trimmed == 'sequences:') {
            inSequences = true
            return
        }
        if (inSequences && !line.startsWith('  ') && !trimmed.isEmpty()) {
            inSequences = false
            return
        }
        if (!inSequences){
            return
        }


        if (trimmed.startsWith('-') && trimmed.endsWith(':')) {
            if (!currentEntry.isEmpty()) {
                enrichedEntries << currentEntry
            }
            currentEntry = ['type': trimmed[1..-2]]
        }else{
            def (key, value) = trimmed.split(':', 2)*.trim()
            currentEntry[key] = value
        }
    }
    if (!currentEntry.isEmpty()) {
        enrichedEntries << currentEntry
    }
    return enrichedEntries
}

def getFastaSequences(fastaData) {
    List<Map> fastaEntries = []
    String currentId = null
    StringBuilder currentSeq = new StringBuilder()

    fastaData.split("\n").each { line ->
        if (line.startsWith(">")) {
            if (currentId) {
                fastaEntries << [id: currentId, sequence: currentSeq.toString()]
            }
            currentId = line[1..-1].trim()  // Remove '>' and trim
            currentSeq = new StringBuilder()
        } else {
            currentSeq.append(line.trim())
        }
    }

    if (currentId) {
        fastaEntries << [id: currentId, sequence: currentSeq.toString()]
    }

    return fastaEntries
}
