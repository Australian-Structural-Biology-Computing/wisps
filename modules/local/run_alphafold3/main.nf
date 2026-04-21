/*
 * Run Alphafold3
 */
process RUN_ALPHAFOLD3 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'
    container "quay.io/nf-core/proteinfold_alphafold3_standard:2.1.0dev"
    input:
    tuple val(meta), path("data/*")
    path "params/*"
    output:
    tuple val(meta), path ("out/*/*_model.cif")       , emit: top_ranked_cif
    tuple val(meta), path ("out/*/*_confidences.json"), emit: jsons
    path "versions.yml"                               , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ALPHAFOLD3 module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 /app/alphafold/run_alphafold.py \\
        --input_dir=data/ \\
        --model_dir=./params \\
        --norun_data_pipeline \\
        --output_dir=out/ \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir out

    for f in data/*; do
        [[ -f "\$f" ]] || continue
        fname=\${f##*/}
        id=\${fname%%.*}
        mkdir -p out/\$id
        touch out/\$id/\${id}_model.cif
        touch out/\$id/\${id}_confidences.json
        touch out/\$id/\${id}_summary_confidences.json
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
