process COLLECT_CONFIDENCE {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(samples)
    path(json_files)
    output:
    tuple val(meta), path ("*confidence_scores_full.csv"), emit: confidence
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env python3

    import os, sys
    import json
    import csv
    import statistics

    # Folder containing the JSON files
    output_csv = "${meta.model}_confidence_scores_full.csv"

    # Collect all JSON files
    json_files = ["${json_files.join('", "')}"]
    samples = ["${samples.join('", "')}"]

    def extract_boltz_results(json_files, samples, output_csv):
        rows = []
        all_keys = ["confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm", "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde", "chains_ptm"]
        writing_keys = ["id", "confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm", "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"]
        writing_keys_extra = set()
        # First pass: read and collect keys
        cntr = 0
        for file in json_files:
            with open(file) as f:
                data = json.load(f)
                flat_data = {"id": samples[cntr]}
                cntr += 1
                for key, value in data.items():
                    if key not in all_keys:
                        continue
                    if isinstance(value, dict):
                        for sub_key, sub_value in value.items():
                            flat_data[f"{key}_{sub_key}"] = sub_value
                            writing_keys_extra.add(f"{key}_{sub_key}")
                    else:
                        flat_data[key] = value
                rows.append(flat_data)

        # Write to CSV

        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=writing_keys + sorted(writing_keys_extra))
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in writing_keys})

    def extract_colabfold_results(json_files, samples, output_csv):
        def parse_colabfold_metrics(json_file):
            with open(json_file, "r") as f:
                data = json.load(f)

            result = {}
            for key in ["iptm", "max_pae", "pae", "ptm", "plddt"]:
                if key not in data:
                    result[key] = None
                    continue

                if key == "plddt" and isinstance(data[key], list):
                    # Compute mean for list
                    result[key] = statistics.mean(data[key]) if data[key] else None
                else:
                    # Take value directly
                    result[key] = data[key]

            return result

        rows = []
        all_keys = ["id", "plddt", "ptm", "iptm", "max_pae"]
        writing_keys = all_keys
        cntr = 0
        for file in json_files:
            flat_data = parse_colabfold_metrics(file)
            flat_data["id"] = samples[cntr]
            cntr += 1
            rows.append(flat_data)

        # Write to CSV
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=writing_keys)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in writing_keys})



    if "${meta.model}".lower() == "boltz":
        extract_boltz_results(json_files, samples, output_csv)
    elif "${meta.model}".lower() == "colabfold":
        extract_colabfold_results(json_files, samples, output_csv)
    elif "${meta.model}".lower() == "alphafold3":
        extract_alphafold3_results(json_files, samples, output_csv)

    print(f"CSV file created at: {output_csv}")

    with open ("versions.yml", "w") as version_file:
        version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """

    stub:
    """
    touch "confidence_scores_full.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        generate_comparison_report.py: \$(python3 --version)
    END_VERSIONS
    """
}
