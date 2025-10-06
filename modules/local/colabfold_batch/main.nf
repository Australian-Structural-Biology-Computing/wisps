process COLABFOLD_BATCH {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    container "nf-core/proteinfold_colabfold:dev"

    input:
    tuple val(meta), path("alignment/*")
    val   colabfold_model_preset
    path  ('params/*')
    path  ('colabfold_db/*')
    path  ('uniref30/*')
    val   numRec

    output:
    tuple val(meta), path ("*relaxed_rank*_model_1_*.pdb")     , emit: pdb
    tuple val(meta), path ("*_coverage.png")          , emit: msa
    tuple val(meta), path ("*_scores_rank_001_*.json")    , emit: top_ranked_scores
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local COLABFOLD_BATCH module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''

    """
    ln -s \$(realpath params/alphafold_params_*/*) params/
    touch params/download_finished.txt

    colabfold_batch \\
        $args \\
        --num-recycle ${numRec} \\
        --data \$PWD \\
        --model-type ${colabfold_model_preset} \\
        alignment \\
        \$PWD
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_batch: \$(pip list | grep "^colabfold" | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ./${meta.id}_relaxed_rank_model_1_00.pdb
    touch ./${meta.id}_relaxed_rank_model_2_00.pdb
    touch ./${meta.id}_relaxed_rank_model_3_00.pdb
    touch ./${meta.id}_coverage.png
    touch ./${meta.id}_scores_rank_001_00.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_batch: \$(pip list | grep "^colabfold" | awk '{print \$2}')
    END_VERSIONS
    """
}
