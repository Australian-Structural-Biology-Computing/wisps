process SPLIT_MSA {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgikit:2.1.0--h3ab6199_0' :
        'biocontainers/mgikit:2.1.0--h3ab6199_0' }"

    input:
    tuple val(meta), path("input_msa/*")
    output:
    tuple val(meta), path ("output_msa/*.csv"), emit: msa_csv
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    msakit split-a3m --input input_msa --output output_msa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msakit: 0.0.1
    END_VERSIONS
    """

    stub:
    """
    mkdir output_msa
    touch "output_msa/A.csv"
    touch "output_msa/B.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msakit: 0.0.1
    END_VERSIONS
    """
}
