process COLLECT_CONFIDENCE {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(samples)
    path("input_json_files/*")
    output:
    tuple val(meta), path ("*confidence_scores_full.csv"), emit: confidence
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def samples_json = groovy.json.JsonOutput.toJson(samples)

    """
    cat << 'EOF' > sample_metadata.json
    ${samples_json}
    EOF

    collect_confidence.py \\
        --input input_json_files \\
        --model "${meta.model}"\\
        --sample-metadata sample_metadata.json \\
        --output "${meta.model}_confidence_scores_full.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        collect_confidence.py: \$(python3 --version)
    END_VERSIONS
    """

    stub:
    """
    touch "confidence_scores_full.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        collect_confidence.py: \$(python3 --version)
    END_VERSIONS
    """
}
