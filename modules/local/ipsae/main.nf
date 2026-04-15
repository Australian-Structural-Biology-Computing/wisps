process IPSAE {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/numpy:2.2.2' :
        'quay.io/biocontainers/numpy:2.2.2' }"

    input:
    tuple val(meta), val(files_names)
    path ("data/")

    output:
    tuple val(meta), path ("data/*[0-9].txt"), emit: txt
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    pae_files=(${files_names.collect { "\"${it[0]}\"" }.join(" ")})
    pdb_files=(${files_names.collect { "\"${it[1]}\"" }.join(" ")})

    [[ \${#pae_files[@]} -eq \${#pdb_files[@]} ]] || { echo "lengths differ"; exit 1; }

    for i in "\${!pae_files[@]}"; do
        ipsae.py data/\${pae_files[i]} \\
                data/\${pdb_files[i]} \\
                ${args}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    pae_files=(${files_names.collect { "\"${it[0]}\"" }.join(" ")})
    pdb_files=(${files_names.collect { "\"${it[1]}\"" }.join(" ")})

    [[ \${#pae_files[@]} -eq \${#pdb_files[@]} ]] || { echo "lengths differ"; exit 1; }

    for i in "\${!pae_files[@]}"; do
        s="\${pdb_files[i]}"
        touch "data/\${s%%_*}_10.pml"
        touch "data/\${s%%_*}_10.txt"
        touch "data/\${s%%_*}_10_byres.txt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
