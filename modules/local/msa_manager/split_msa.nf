process SPLIT_MSA {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(msa)
    output:
    tuple val(meta), path ("output_msa/*.csv"), emit: msa_csv
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    files=(${msa.collect { "\"${it}\"" }.join(" ")})

    for f in "\${files[@]}"; do
        meta_id=\$(basename "\$f" ".a3m")
        msa_manager.py \$f -o output_msa --meta_id "\${meta_id}"
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir output_msa
    touch "output_msa/A.csv"
    touch "output_msa/B.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
