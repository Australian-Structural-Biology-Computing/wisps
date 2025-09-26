process PREPARE_INTERACTIONS {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(fasta)
    tuple val(meta), val(a3m)
    tuple val(meta), val(types)
    path(files)
    output:
    tuple val(meta), path ("*.fasta"), emit: fasta
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env python3
    import os, sys
    import string
    
    #single_letters = list(string.ascii_uppercase)
    #two_letter_combinations = [''.join(p) for p in product(string.ascii_uppercase, repeat=2)]
    #all_combinations = single_letters + two_letter_combinations
    all_combinations = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [str(x) for x in range(0, 10)]
    fasta_files = ["${fasta.join('", "')}"]
    a3m_files = ["${a3m.join('", "')}"]
    seq_types = ["${types.join('", "')}"]
    for seq_type in seq_types:
        if seq_type not in {"protein", "dna", "rna", "smiles", "ccd"}:
            print(f'seqeuence type should be ["protein", "dna", "rna", "smiles", "ccd"], found {seq_type}!')
            exit (1)
    output_file = "${meta.id}_boltz_interaction_input.fasta"
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
