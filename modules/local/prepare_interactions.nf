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
        tuple val(meta), path ("*_boltz_interaction_input.fasta"), emit: fasta
        path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env python3
    import os, sys
    import string
    ids_batch = ${ids.collect { "\"${it}\"" }.toString()}
    fasta_files_batch = ${fasta.collect { inner -> inner.collect { "\"${it}\"" }}.toString()}
    a3m_files_batch = ${a3m.collect { inner -> inner.collect { "\"${it}\"" }}.toString()}
    seq_types_batch = ${types.collect { inner -> inner.collect { "\"${it}\"" }}.toString()}
    
    all_combinations = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [str(x) for x in range(0, 10)]
    batch_itr = 0
    for batch_itr in range(len(fasta_files_batch)):
        a3m_files = a3m_files_batch[batch_itr]
        seq_types = seq_types_batch[batch_itr]
        fasta_files = fasta_files_batch[batch_itr]
        
        for seq_type in seq_types:
            if seq_type not in {"protein", "dna", "rna", "smiles", "ccd"}:
                print(f'seqeuence type should be ["protein", "dna", "rna", "smiles", "ccd"], found {seq_type}!')
                exit (1)
        output_file = ids_batch[batch_itr] + "_boltz_interaction_input.fasta"
        counter = 0
        with open(output_file, "w") as outfile:
            for seq_itr in range(len(seq_types)):
                fasta = fasta_files[seq_itr]
                a3m   = a3m_files[seq_itr]
                seq_type  = seq_types[seq_itr]

                with open(fasta, "r") as f:
                    lines = f.readlines()

                if not lines:
                    continue  # Skip empty FASTA files

                #header = lines[0].strip()
                lines[-1] = lines[-1].strip() + "\\n"
                body = lines[1:]
                header = f">{all_combinations[counter]}|{seq_type}"
                counter += 1
                if seq_type == "protein":
                    header += f"|{os.path.basename(a3m)}\\n"
                else:
                    header += "\\n"
                
                outfile.write(header)
                outfile.writelines(body)

    with open ("versions.yml", "w") as version_file:
	    version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """

    stub:
    """
    touch "${meta.id}_boltz_interaction_input.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        generate_comparison_report.py: \$(python3 --version)
    END_VERSIONS
    """
}
