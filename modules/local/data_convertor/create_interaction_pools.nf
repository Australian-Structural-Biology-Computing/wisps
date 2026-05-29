process CREATE_INTERACTION_POOLS {
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
    val(pool_mode)
    val(pool_max_total_length)

    output:
    path ("interactions.fasta"), emit: interactions
    path ("interaction_mapping.tsv"), emit: interaction_mapping
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import sys

    seqs = [${"\"" + seqs.join('", "') + "\""}]
    ids = [${"\"" + ids.join('", "') + "\""}]
    types = [${"\"" + types.join('", "') + "\""}]
    groups = [${"\"" + groups.join('", "') + "\""}]
    pool_mode = [${"\"" + pool_mode.join('", "') + "\""}]

    pool_max_total_length = int(${pool_max_total_length})

    if pool_max_total_length <= 0:
        raise ValueError(f"pool_max_total_length must be > 0, got {pool_max_total_length}")

    n = len(seqs)
    interactions_cntr = 0
    supported = {"dna", "rna", "ccd", "smiles"}

    def format_non_protein_entry(entity_type, entity_seq):
        parsed_type = entity_type.lower()
        parsed_seq = entity_seq

        # If sequence is already in colabfold format, keep it and normalize type.
        if "|" in entity_seq:
            parts = entity_seq.split("|")
            if len(parts) > 1 and parts[0].lower() in supported:
                parsed_type = parts[0].lower()
                parsed_seq = "|".join(parts[1:])

        if parsed_type == "smiles":
            # ColabFold expects aromatic bonds to be written with ';' instead of ':'.
            if "|" in parsed_seq:
                smi_parts = parsed_seq.split("|")
                smi_parts[0] = smi_parts[0].replace(":", ";")
                parsed_seq = "|".join(smi_parts)
            else:
                parsed_seq = parsed_seq.replace(":", ";")

        return f"{parsed_type}|{parsed_seq}"

    def chain_count(entity_type, entity_seq):
        if entity_type.lower() == "protein":
            return len([chain for chain in entity_seq.split(":") if chain])
        return 1

    def sequence_length(entity_type, entity_seq):
        parsed_type = entity_type.lower()
        parsed_seq = entity_seq

        # Handle already-prefixed non-protein entries: type|seq(|copies)
        if "|" in entity_seq:
            parts = entity_seq.split("|")
            if len(parts) > 1 and parts[0].lower() in supported:
                parsed_type = parts[0].lower()
                parsed_seq = "|".join(parts[1:])

        if parsed_type == "protein":
            return sum(len(chain) for chain in parsed_seq.split(":") if chain)

        # For non-protein entities, use payload length and optional copies.
        non_protein_parts = parsed_seq.split("|")
        payload = non_protein_parts[0]
        copies = 1
        if len(non_protein_parts) > 1 and non_protein_parts[-1].isdigit():
            copies = int(non_protein_parts[-1])
        return len(payload) * copies

    def entity_to_fasta(entity_type, entity_seq):
        if entity_type.lower() == "protein":
            return entity_seq
        return format_non_protein_entry(entity_type, entity_seq)

    def int_id_to_str_id(i):
        if i < 0:
            raise ValueError(f"int_id_to_str_id: Only positive integers allowed, got {i}")
        output = []
        while i >= 0:
            output.append(chr(i % 26 + ord("A")))
            i = i // 26 - 1
        return "".join(output)

    def best_fit_decreasing_pool(indices, lengths, capacity):
        # Bin-packing heuristic: keeps pools close to capacity without going over.
        items = sorted(indices, key=lambda idx: lengths[idx], reverse=True)
        bins = []

        for idx in items:
            item_len = lengths[idx]
            if item_len > capacity:
                raise ValueError(
                    f"Sequence '{ids[idx]}' length {item_len} exceeds available pool capacity {capacity}."
                )

            best_bin_pos = -1
            best_remaining_after = None
            for b_idx, b in enumerate(bins):
                if b["remaining"] >= item_len:
                    rem_after = b["remaining"] - item_len
                    if best_remaining_after is None or rem_after < best_remaining_after:
                        best_remaining_after = rem_after
                        best_bin_pos = b_idx

            if best_bin_pos == -1:
                bins.append({"items": [idx], "remaining": capacity - item_len})
            else:
                bins[best_bin_pos]["items"].append(idx)
                bins[best_bin_pos]["remaining"] -= item_len

        return [b["items"] for b in bins]

    pool_pairs = []
    for mode_spec in pool_mode:
        parts = mode_spec.split("-")
        if len(parts) != 2 or any(not p.strip() for p in parts):
            raise ValueError(
                f"Invalid pooled mode entry '{mode_spec}'. Expected '<group1>-<group2>'."
            )
        pool_pairs.append((parts[0].strip(), parts[1].strip()))

    seq_lengths = [sequence_length(types[i], seqs[i]) for i in range(n)]

    with open("interactions.fasta", "w") as target_interactions, open("interaction_mapping.tsv", "w") as mapping_file:
        mapping_file.write("interaction_id\\tint_chain_map\\tstr_chain_map\\treport_id\\tleft_source_id\\tright_source_id\\n")

        for pool_group_1, pool_group_2 in pool_pairs:
            left_indices = [i for i in range(n) if groups[i] == pool_group_1]
            right_indices = [i for i in range(n) if groups[i] == pool_group_2]

            if not left_indices:
                raise ValueError(f"No entries found in pool group '{pool_group_1}'.")
            if not right_indices:
                raise ValueError(f"No entries found in pool group '{pool_group_2}'.")

            for left_i in left_indices:
                left_len = seq_lengths[left_i]
                pool_capacity = pool_max_total_length - left_len
                if pool_capacity <= 0:
                    raise ValueError(
                        f"Left sequence '{ids[left_i]}' length {left_len} is >= pool_max_total_length {pool_max_total_length}."
                    )

                # If groups are identical, do not self-pool with the same row.
                current_right_indices = [idx for idx in right_indices if idx != left_i]
                if not current_right_indices:
                    raise ValueError(
                        f"No valid group2 entries available for left entry '{ids[left_i]}' after self-exclusion."
                    )

                pools = best_fit_decreasing_pool(current_right_indices, seq_lengths, pool_capacity)

                for pool_idx, pooled_right in enumerate(pools, start=1):
                    # Keep only pools that include at least one protein partner overall.
                    has_left_protein = types[left_i].lower() == "protein"
                    has_any_right_protein = any(types[r].lower() == "protein" for r in pooled_right)
                    if not (has_left_protein or has_any_right_protein):
                        continue

                    left_chain_count = chain_count(types[left_i], seqs[left_i])
                    right_chain_counts = [chain_count(types[r], seqs[r]) for r in pooled_right]

                    left_chain_indices = ",".join(str(idx) for idx in range(0, left_chain_count))
                    left_chain_str_ids = ",".join(int_id_to_str_id(idx) for idx in range(0, left_chain_count))

                    interaction_id = f"{pool_group_1}-{pool_group_2}__{ids[left_i]}-pool{pool_idx:03d}"

                    entities = [entity_to_fasta(types[left_i], seqs[left_i])]
                    entities.extend(entity_to_fasta(types[r], seqs[r]) for r in pooled_right)
                    target_interactions.write(f">{interaction_id}\\n{':'.join(entities)}\\n")
                    interactions_cntr += 1

                    # Emit one mapping row per original cross-group pair in this pool.
                    right_start = left_chain_count
                    for ridx, right_i in enumerate(pooled_right):
                        right_count = right_chain_counts[ridx]
                        right_indices_pair = ",".join(
                            str(idx) for idx in range(right_start, right_start + right_count)
                        )
                        right_str_ids_pair = ",".join(
                            int_id_to_str_id(idx) for idx in range(right_start, right_start + right_count)
                        )
                        report_id = f"{ids[left_i]}-{ids[right_i]}"
                        mapping_file.write(
                            f"{interaction_id}\\t"
                            f"{left_chain_indices}:{right_indices_pair}\\t"
                            f"{left_chain_str_ids}:{right_str_ids_pair}\\t"
                            f"{report_id}\\t{ids[left_i]}\\t{ids[right_i]}\\n"
                        )
                        right_start += right_count

    if interactions_cntr == 0:
        print("No interaction pools found!!")
        exit(1)

    print(f"{interactions_cntr} Interaction pools found!!")
    with open ("versions.yml", "w") as version_file:
        version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """

    stub:
    """
    touch interactions.fasta
    touch interaction_mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
