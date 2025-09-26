process IPSAE {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/numpy:2.2.2' :
        'quay.io/biocontainers/numpy:2.2.2' }"

    input:
    tuple val(meta), path(pae), path(pdb)
    output:
    tuple val(meta), path ("*.pml"), emit: pml
    tuple val(meta), path ("*_byres.txt"), emit: byres
    tuple val(meta), path ("*.txt"), emit: txt
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    ipsae.py ${pae} \\
        ${pdb} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.pml"
    touch "${meta.id}.txt"
    touch "${meta.id}.byres.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
