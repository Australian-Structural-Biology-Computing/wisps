process MMSEQS_COLABFOLDSEARCH {
    tag "$meta.id"
    label 'process_high_memory'

    container "ghcr.io/tlitfin/wisps-colabfold-search:1.1"

    input:
    tuple val(meta), path(fasta)
    path ('db/*')
    path ('db/*')

    output:
    tuple val(meta), path("**.a3m"), emit: a3m
    tuple val(meta), path("**.json"), emit: json, optional: true
    tuple val(meta), path("**.yaml"), emit: yaml, optional: true
    tuple val(meta), path("**.csv"), emit: msa_csv, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local MMSEQS_COLABFOLDSEARCH module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''

    """
    colabfold_search \\
        $args \\
        --threads $task.cpus \\
        ${fasta} \\
        ./db \\
        "result/"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_search: \$(pip list | grep "^colabfold" | awk '{print \$2}')
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    """
    mkdir results

    if [[ "$fasta" == *.csv ]]; then
        line_num=0
        while IFS=, read -r firstcol rest; do
            line_num=\$((line_num + 1))
            # Skip the first line (header)
            if [ "\$line_num" -eq 1 ]; then
            continue
            fi
            # Skip empty first column
            [ -z "\$firstcol" ] && continue

            # Clean filename (letters, numbers, dash, underscore)
            safe_name=\$(echo "\$firstcol" | tr -cd '[:alnum:]_-')

            touch "results/\${safe_name}.a3m"
            touch "results/\${safe_name}.json"
        done < "$fasta"
    else
        touch results/${meta.id}.a3m
        touch results/${meta.id}.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_search: \$(conda run -n colabfold pip list | grep "^colabfold" | awk '{print \$2}')
        alphafold_colabfold: \$(conda run -n colabfold pip list | grep "^alphafold-colabfold" | awk '{print \$2}')
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """
}
