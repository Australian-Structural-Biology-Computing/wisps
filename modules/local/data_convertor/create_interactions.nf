process CREATE_INTERACTIONS {
    tag   "all_run"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    val(ids)
    val(types)
    val(groups)
    val(seqs)
    val(interaction_mode)
    val(interaction_threshold)

    output:
    path ("interactions.fasta"), emit: interactions
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env python3
    import os, sys
    interaction_mode = [${"\"" + interaction_mode.join('", "') + "\""}]
    seqs = [${"\"" + seqs.join('", "') + "\""}]
    ids = [${"\"" + ids.join('", "') + "\""}]
    types = [${"\"" + types.join('", "') + "\""}]
    groups = [${"\"" + groups.join('", "') + "\""}]
    done_singles = set()
    n = len(seqs)
    interactions_cntr = 0
    group_ls = {part for s in groups for part in s.split("-")}
    with open("interactions.fasta", "w") as target_interactions:
        for i in range(n):
            for j in range(i, n):
                if interaction_mode[0] != "all-all" and groups[i] not in group_ls:
                    continue

                if (interaction_mode[0] != "all-all" and "-".join(sorted([groups[i], groups[j]])) in interaction_mode) or (interaction_mode[0] == "all-all" and (${interaction_threshold} == 0 or abs(i-j) <= ${interaction_threshold})):
                    interactions_cntr += 1
                    if types[i] == "protein" and types[j] == "protein":
                        target_interactions.write(f">{ids[i]}-{ids[j]}\\n{seqs[i]}:{seqs[j]}\\n")
                    elif types[i] == "protein":
                        if ids[i] not in done_singles:
                            target_interactions.write(f">{ids[i]}\\n{seqs[i]}\\n")
                            done_singles.add(ids[i])
                    elif types[j] == "protein":
                        if ids[j] not in done_singles:
                            target_interactions.write(f">{ids[j]}\\n{seqs[j]}\\n")
                            done_singles.add(ids[j])
    if interactions_cntr == 0:
        print("No Interactions found!!")
        exit(1)
    print(f"{interactions_cntr} Interactions found!!")
    with open ("versions.yml", "w") as version_file:
        version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """

    stub:
    """
    touch interactions.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
