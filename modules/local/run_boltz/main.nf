/*
 * Run Boltz
 */
process RUN_BOLTZ {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    container "nf-core/proteinfold_boltz:dev"

    input:
    tuple val(meta), path("fasta/*")
    path(alignments)
    path ('boltz1_conf.ckpt')
    path ('ccd.pkl')
    path ('boltz2_aff.ckpt')
    path ('boltz2_conf.ckpt')
    path ('mols')

    output:
    tuple val(meta), path ("boltz_results_*/processed/msa/*.npz")               , optional: true, emit: msa
    tuple val(meta), path ("boltz_results_*/processed/structures/*.npz")        , emit: structures
    tuple val(meta), path ("boltz_results_*/predictions/*/confidence*.json")    , emit: confidence
    tuple val(meta), path ("boltz_results_*/predictions/*/*.pdb")               , optional: true, emit: pdb
    tuple val(meta), path ("boltz_results_*/predictions/*/*model_0.cif")        , optional: true, emit: cif
    tuple val(meta), path ("boltz_results_*/predictions/*/plddt_*model_0.npz")  , emit: plddt
    tuple val(meta), path ("boltz_results_*/predictions/*/pae_*model_0.npz")    , optional: true, emit: pae
    
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_BOLTZ module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def version = "2.0.3"
    def args = task.ext.args ?: ''

    """
    export NUMBA_CACHE_DIR=/tmp
    export HOME=/tmp

    boltz predict fasta ${args} --cache ./
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        boltz: $version
    END_VERSIONS
    """

    stub:
    def version = "2.0.3"
    """
    mkdir -p boltz_results_${meta.id}/processed/msa/
    mkdir -p boltz_results_${meta.id}/processed/structures/
    mkdir -p boltz_results_${meta.id}/predictions/${meta.id}/

    touch boltz_results_${meta.id}/processed/msa/${meta.id}.npz
    touch boltz_results_${meta.id}/processed/structures/${meta.id}.npz
    touch boltz_results_${meta.id}/predictions/${meta.id}/confidence_${meta.id}.json
    touch boltz_results_${meta.id}/predictions/${meta.id}/${meta.id}.pdb
    touch boltz_results_${meta.id}/predictions/${meta.id}/${meta.id}.cif
    touch boltz_results_${meta.id}/predictions/${meta.id}/plddt_${meta.id}_model_0.npz
    touch boltz_results_${meta.id}/predictions/${meta.id}/pae_${meta.id}_model_0.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        boltz: $version
    END_VERSIONS
    """
}
