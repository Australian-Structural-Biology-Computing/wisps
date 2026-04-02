process PREPARE_INTERACTIONS {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
        tuple val(meta), val(ids)
        val(fasta)
        val(a3m)
        val(types)
        path(files)

    output:
        tuple val(meta), path ("*_boltz_interaction_input.yaml"), emit: yaml
        path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env python3
    import os, sys
    import string
    import json

    ids_batch = ${ids.collect { "\"${it}\"" }.toString()}
    fasta_files_batch = ${fasta.collect { inner -> inner.collect { "\"${it}\"" }}.toString()}
    a3m_files_batch = ${a3m.collect { inner -> inner.collect { "\"${it}\"" }}.toString()}
    seq_types_batch = ${types.collect { inner -> inner.collect { "\"${it}\"" }}.toString()}

    all_combinations = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [str(x) for x in range(0, 10)]

    def read_fasta_sequence(fasta_path):
        with open(fasta_path, "r") as f:
            lines = [line.strip() for line in f.readlines()]

        if not lines:
            return ""

        sequence_lines = [line for line in lines if line and not line.startswith(">")]
        return "".join(sequence_lines)

    for batch_itr in range(len(fasta_files_batch)):
        a3m_files = a3m_files_batch[batch_itr]
        seq_types = seq_types_batch[batch_itr]
        fasta_files = fasta_files_batch[batch_itr]

        for seq_type in seq_types:
            if seq_type not in {"protein", "dna", "rna", "smiles", "ccd"}:
                print(f'seqeuence type should be ["protein", "dna", "rna", "smiles", "ccd"], found {seq_type}!')
                exit(1)

        output_file = ids_batch[batch_itr] + "_boltz_interaction_input.yaml"
        counter = 0
        yaml_lines = ["version: 1", "sequences:"]
        ligand_chain_ids = []

        for seq_itr in range(len(seq_types)):
            fasta_path = fasta_files[seq_itr]
            a3m_path = a3m_files[seq_itr]
            seq_type = seq_types[seq_itr]

            sequence = read_fasta_sequence(fasta_path)
            if not sequence:
                continue

            chain_id = all_combinations[counter]
            counter += 1

            if seq_type in {"smiles", "ccd"}:
                yaml_lines.append("  - ligand:")
                yaml_lines.append(f"      id: {chain_id}")
                yaml_lines.append(f"      {seq_type}: {json.dumps(sequence)}")
                ligand_chain_ids.append(chain_id)
            else:
                yaml_lines.append(f"  - {seq_type}:")
                yaml_lines.append(f"      id: {chain_id}")
                yaml_lines.append(f"      sequence: {json.dumps(sequence)}")
                if seq_type == "protein" and a3m_path:
                    yaml_lines.append(f"      msa: {json.dumps(os.path.basename(a3m_path))}")

        if len(ligand_chain_ids) == 1:
            yaml_lines.append("properties:")
            yaml_lines.append("  - affinity:")
            yaml_lines.append(f"      binder: {ligand_chain_ids[0]}")

        with open(output_file, "w") as outfile:
            outfile.write("\\n".join(yaml_lines) + "\\n")

    with open("versions.yml", "w") as version_file:
        version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """

    stub:
    """
    #!/usr/bin/env python3
    from pathlib import Path
    import sys

    ids_batch = ${ids.collect { "\"${it}\"" }.toString()}
    for batch_itr in range(len(ids_batch)):
        Path(ids_batch[batch_itr] + "_boltz_interaction_input.yaml").touch()

    with open("versions.yml", "w") as version_file:
        version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """
}
