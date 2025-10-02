/*
 * Run Alphafold3
 */
process RUN_ALPHAFOLD3 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'
    container "<path to alphafold3.sif>"
    input:
    tuple val(meta), path(json)
    path "params/*"
    output:
    tuple val(meta), path ("publish/*alphafold3.cif")       , emit: top_ranked_cif
    tuple val(meta), path ("publish/*ranked_*.cif")         , emit: cif
    tuple val(meta), path ("${meta.id}_plddt.tsv")          , emit: multiqc
    tuple val(meta), path ("${meta.id}_alphafold3_msa.tsv") , emit: msa
    tuple val(meta), path ("${meta.id}_0_pae.tsv")          , emit: pae
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
        --json_path=${json} \\
        --model_dir=./params \\
        --norun_data_pipeline \\
        --output_dir=\$PWD \\
        $args
    ## Rename the top ranked model
    if [ ! -d publish ]; then
        mkdir -p publish
    fi
    ## Move the rest of the models and rename them according to their rank
    name=\$(jq -r '.name' ${json})
    cp -n "\${name}/\${name}_model.cif" "publish/${prefix}_alphafold3.cif"
    # Sort the rows by ranking_score in descending order
    sorted_csv=\$(head -n 1 "\${name}/\${name}_ranking_scores.csv"; tail -n +2 "\${name}/\${name}_ranking_scores.csv" | sort -t, -k3 -nr)
    rank=0
    touch publish/combined_plddt_mqc.tsv
    # Generate files with rank tag
    echo "\$sorted_csv" | tail -n +2 | while IFS=',' read -r seed sample ranking_score; do
    cp -n "\${name}/seed-\${seed}_sample-\${sample}/\${name}_seed-\${seed}_sample-\${sample}_model.cif" "publish/seed_\${seed}_sample_\${sample}_ranked_\${rank}.cif"
    rank=\$((rank + 1))
    done
    extract_metrics.py --name ${prefix} \\
        --jsons ${meta.id}/${meta.id}_data.json ${meta.id}/${meta.id}_summary_confidences.json ${meta.id}/${meta.id}_confidences.json \\
        --structs publish/*ranked_*.cif
    mv "${meta.id}_msa.tsv" "${meta.id}_alphafold3_msa.tsv"
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
