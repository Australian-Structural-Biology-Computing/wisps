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
    tuple val(meta), path ("publish/*alphafold3.cif")       , emit: top_ranked_cif
    tuple val(meta), path ("publish/*ranked_*.cif")         , emit: cif
    tuple val(meta), path ("${meta.id}/${meta.id}_confidences.json")   , emit: raw_pae
    path "versions.yml"                                     , emit: versions
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
        --input_dir=${data} \\
        --model_dir=./params \\
        --norun_data_pipeline \\
        --output_dir=\$PWD \\
        $args
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir publish
    touch publish/${prefix}_alphafold3.cif
    touch publish/${prefix}_ranked_1.cif
    touch publish/${prefix}_ranked_2.cif
    touch publish/${prefix}_ranked_3.cif
    touch publish/${prefix}_ranked_4.cif
    touch publish/${prefix}_ranked_5.cif
    touch ${prefix}_plddt.tsv
    touch ${prefix}_alphafold3_msa.tsv
    touch ${prefix}_0_pae.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
