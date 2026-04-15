#!/usr/bin/env python3

import json
import csv
import statistics
import re
from pathlib import Path
from typing import Any, Dict, Iterable
import argparse

def flatten_deep(
    data: Dict[Any, Any],
    *,
    sep: str = "_",
    index_lists: bool = False,   # set True to flatten lists/tuples by index
    keep_empty_dicts: bool = False
) -> Dict[str, Any]:
    """
    Recursively flattens nested dictionaries by joining keys with `sep`.
    - Dict keys at any depth are converted to strings in the output keys.
    - If `index_lists` is True, lists/tuples are traversed and indexed.
    - If `keep_empty_dicts` is True, empty dicts become key -> {} (else skipped).
    """
    flat: Dict[str, Any] = {}

    def _walk(prefix: str, obj: Any) -> None:
        # dict -> dive deeper
        if isinstance(obj, dict):
            if not obj:
                if keep_empty_dicts and prefix:
                    flat[prefix] = {}
                return
            for k, v in obj.items():
                k_str = str(k)
                new_key = f"{prefix}{sep}{k_str}" if prefix else k_str
                _walk(new_key, v)
            return

        # list/tuple handling (optional)
        if index_lists and isinstance(obj, (list, tuple)):
            if not obj:
                # represent empty list if desired; default: skip
                return
            for i, item in enumerate(obj):
                new_key = f"{prefix}{sep}{i}" if prefix else str(i)
                _walk(new_key, item)
            return

        # leaf value -> record it
        if prefix == "":
            # top-level non-dict without a key path; use a default name if needed
            flat["_"] = obj
        else:
            flat[prefix] = obj

    _walk("", data)
    return flat

def parse_boltz_json_file(json_file):
    name = Path(json_file).name
    if name.startswith("confidence_") and "_boltz_" in name:
        return "confidence", name.split("_boltz_")[0].split("confidence_")[1]
    if name.startswith("affinity_") and "_boltz_" in name:
        return "affinity", name.split("_boltz_")[0].split("affinity_")[1]

    m = re.match(r"confidence_(.+?)_model_\d+\.json$", name)
    if m:
        return "confidence", m.group(1)
    m = re.match(r"affinity_(.+?)\.json$", name)
    if m:
        return "affinity", m.group(1)
    return None, None

def build_sample_meta_by_id(sample_metadata):
    sample_meta_by_id = {}
    if not sample_metadata:
        return sample_meta_by_id

    for item in sample_metadata:
        if isinstance(item, dict) and "id" in item:
            sample_id = str(item["id"])
            sample_meta_by_id.setdefault(sample_id, []).append({
                "id": sample_id,
                "report_id": item.get("report_id", sample_id),
                "int_chain_map": item.get("int_chain_map", ""),
                "str_chain_map": item.get("str_chain_map", ""),
                "left_source_id": item.get("left_source_id", ""),
                "right_source_id": item.get("right_source_id", ""),
            })
    return sample_meta_by_id

def get_meta_entries(sample_id, sample_meta_by_id):
    entries = sample_meta_by_id.get(sample_id, [])
    if entries:
        return entries
    return [{
        "id": sample_id,
        "report_id": sample_id,
        "int_chain_map": "",
        "str_chain_map": "",
        "left_source_id": "",
        "right_source_id": "",
    }]

def apply_mapping_fields(row, meta):
    row["id"] = meta.get("report_id", row.get("id", ""))
    row["model_input_id"] = meta.get("id", "")
    row["int_chain_map"] = meta.get("int_chain_map", "")
    row["str_chain_map"] = meta.get("str_chain_map", "")
    row["left_source_id"] = meta.get("left_source_id", "")
    row["right_source_id"] = meta.get("right_source_id", "")

def is_empty_value(value):
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    return False

def prune_empty_columns(rows, fieldnames):
    retained = []
    for field in fieldnames:
        if any(not is_empty_value(row.get(field)) for row in rows):
            retained.append(field)
    return retained

def format_row_for_multiqc_csv(row, fieldnames):
    """
    Emit deterministic table rows for MultiQC:
    - no empty cells (use 'NA')
    - avoid commas in chain-map display fields
    """
    formatted = {}
    for field in fieldnames:
        value = row.get(field, "")
        if is_empty_value(value):
            value = "NA"
        if field in {"int_chain_map", "str_chain_map"} and isinstance(value, str):
            value = value.replace(",", ";")
        formatted[field] = value
    return formatted

def parse_int_chain_map(int_chain_map):
    if not int_chain_map or ":" not in int_chain_map:
        return [], []

    left, right = int_chain_map.split(":", 1)
    left_ids = [int(x) for x in left.split(",") if x.strip() != ""]
    right_ids = [int(x) for x in right.split(",") if x.strip() != ""]
    return left_ids, right_ids

def parse_str_chain_map(str_chain_map):
    if not str_chain_map or ":" not in str_chain_map:
        return [], []

    left, right = str_chain_map.split(":", 1)
    left_ids = [x.strip() for x in left.split(",") if x.strip() != ""]
    right_ids = [x.strip() for x in right.split(",") if x.strip() != ""]
    return left_ids, right_ids

def compute_max_cross_group_metric(int_chain_map, value_getter):
    left_ids, right_ids = parse_int_chain_map(int_chain_map)
    if not left_ids or not right_ids:
        return None

    values = []
    for li in left_ids:
        for rj in right_ids:
            v_lr = value_getter(li, rj)
            v_rl = value_getter(rj, li)

            if isinstance(v_lr, (int, float)):
                values.append(v_lr)
            if isinstance(v_rl, (int, float)):
                values.append(v_rl)

    return max(values) if values else None

def compute_max_cross_group_metric_from_groups(left_ids, right_ids, value_getter):
    if not left_ids or not right_ids:
        return None

    values = []
    for left_id in left_ids:
        for right_id in right_ids:
            v_lr = value_getter(left_id, right_id)
            v_rl = value_getter(right_id, left_id)

            if isinstance(v_lr, (int, float)):
                values.append(v_lr)
            if isinstance(v_rl, (int, float)):
                values.append(v_rl)

    return max(values) if values else None

def compute_overall_max_pair_chains_iptm(pair_chains_iptm):
    if not isinstance(pair_chains_iptm, dict):
        return None

    values = []
    for left_id, right_map in pair_chains_iptm.items():
        if not isinstance(right_map, dict):
            continue
        for right_id, value in right_map.items():
            if str(left_id).strip() == str(right_id).strip():
                continue
            if isinstance(value, (int, float)):
                values.append(value)

    return max(values) if values else None

def compute_overall_max_chain_pair_iptm(chain_pair_iptm):
    if not isinstance(chain_pair_iptm, list):
        return None

    values = []
    for i, row in enumerate(chain_pair_iptm):
        if not isinstance(row, list):
            continue
        for j, value in enumerate(row):
            if i == j:
                continue
            if isinstance(value, (int, float)):
                values.append(value)

    return max(values) if values else None

def split_pair_key(pair_key):
    key = str(pair_key).strip()
    if "-" in key:
        parts = key.split("-", 1)
    elif "_" in key:
        parts = key.split("_", 1)
    else:
        return None, None

    if len(parts) != 2:
        return None, None
    return parts[0].strip(), parts[1].strip()

def compute_overall_max_pairwise_iptm(pairwise_iptm):
    if not isinstance(pairwise_iptm, dict):
        return None

    values = []
    for pair_key, value in pairwise_iptm.items():
        left_id, right_id = split_pair_key(pair_key)
        if left_id is not None and right_id is not None and left_id == right_id:
            continue
        if isinstance(value, (int, float)):
            values.append(value)

    return max(values) if values else None

def compute_max_cross_group_pair_chains_iptm(pair_chains_iptm, int_chain_map):
    if not isinstance(pair_chains_iptm, dict):
        return None
    overall_max = compute_overall_max_pair_chains_iptm(pair_chains_iptm)

    def get_dict_value(i, j):
        return pair_chains_iptm.get(str(i), {}).get(str(j))

    cross_group_max = compute_max_cross_group_metric(int_chain_map, get_dict_value)
    return cross_group_max if cross_group_max is not None else overall_max

def compute_max_cross_group_chain_pair_iptm(chain_pair_iptm, int_chain_map):
    if not isinstance(chain_pair_iptm, list):
        return None
    overall_max = compute_overall_max_chain_pair_iptm(chain_pair_iptm)

    def get_list_value(i, j):
        if 0 <= i < len(chain_pair_iptm) and isinstance(chain_pair_iptm[i], list):
            if 0 <= j < len(chain_pair_iptm[i]):
                return chain_pair_iptm[i][j]
        return None

    cross_group_max = compute_max_cross_group_metric(int_chain_map, get_list_value)
    return cross_group_max if cross_group_max is not None else overall_max

def compute_max_cross_group_pairwise_iptm(pairwise_iptm, str_chain_map):
    if not isinstance(pairwise_iptm, dict):
        return None
    overall_max = compute_overall_max_pairwise_iptm(pairwise_iptm)

    left_ids, right_ids = parse_str_chain_map(str_chain_map)

    def get_pairwise_value(i, j):
        return pairwise_iptm.get(f"{i}-{j}")

    cross_group_max = compute_max_cross_group_metric_from_groups(left_ids, right_ids, get_pairwise_value)
    return cross_group_max if cross_group_max is not None else overall_max

def extract_boltz_results(json_files, output_csv, sample_meta_by_id):
    rows_by_id = {}
    all_keys = ["confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm", "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde", "chains_ptm"]
    writing_keys = ["id", "model_input_id", "left_source_id", "right_source_id", "iptm", "ptm", "int_chain_map", "str_chain_map", "max_cross_group_pair_chains_iptm", "confidence_score", "ligand_iptm", "protein_iptm", "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"]
    writing_keys_extra = set()

    for file in json_files:
        record_type, sample_id = parse_boltz_json_file(file)
        if sample_id is None:
            continue

        if sample_id not in rows_by_id:
            rows_by_id[sample_id] = {"id": sample_id}

        with open(file) as f:
            data = json.load(f)
            flat_data = rows_by_id[sample_id]

            if record_type == "confidence":
                flat_data["_pair_chains_iptm"] = data.get("pair_chains_iptm")
                for key, value in data.items():
                    if key not in all_keys:
                        continue
                    if isinstance(value, dict):
                        for sub_key, sub_value in value.items():
                            key_name = f"{key}_{sub_key}"
                            flat_data[key_name] = sub_value
                            writing_keys_extra.add(key_name)
                    else:
                        flat_data[key] = value

            elif record_type == "affinity":
                for key, value in data.items():
                    if isinstance(value, dict) or isinstance(value, list):
                        sub_data = flatten_deep(value, index_lists=True)
                        for sub_key, sub_value in sub_data.items():
                            key_name = f"{key}_{sub_key}"
                            if not key_name.startswith("affinity_"):
                                key_name = f"affinity_{key_name}"
                            flat_data[key_name] = sub_value
                            writing_keys_extra.add(key_name)
                    else:
                        key_name = key if key.startswith("affinity_") else f"affinity_{key}"
                        flat_data[key_name] = value
                        writing_keys_extra.add(key_name)

    # Prepare output rows and prune columns with no data
    final_writing_vals = writing_keys + sorted(writing_keys_extra - set(writing_keys))
    output_rows = []

    for sample_id in sorted(rows_by_id.keys()):
        base_row = rows_by_id[sample_id]
        for meta in get_meta_entries(sample_id, sample_meta_by_id):
            row = dict(base_row)
            apply_mapping_fields(row, meta)
            row["max_cross_group_pair_chains_iptm"] = compute_max_cross_group_pair_chains_iptm(
                row.get("_pair_chains_iptm"),
                row.get("int_chain_map", "")
            )
            output_rows.append({key: row.get(key, "") for key in final_writing_vals})

    final_writing_vals = prune_empty_columns(
        output_rows,
        final_writing_vals
    )

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=final_writing_vals)
        writer.writeheader()
        for row in output_rows:
            writer.writerow(format_row_for_multiqc_csv(row, final_writing_vals))

def extract_colabfold_results(json_files, samples, output_csv, sample_meta_by_id):
    def parse_colabfold_metrics(json_file):
        with open(json_file, "r") as f:
            data = json.load(f)

        result = {}
        for key in ["iptm", "pae", "ptm", "plddt", "pairwise_iptm"]:
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
    all_keys = ["id", "model_input_id", "left_source_id", "right_source_id", "ptm", "iptm", "int_chain_map", "str_chain_map", "max_cross_group_pairwise_iptm", "plddt"]
    writing_keys = all_keys
    cntr = 0
    for file in json_files:
        sample_id = samples[cntr]
        base_data = parse_colabfold_metrics(file)
        base_data["id"] = sample_id
        cntr += 1
        for meta in get_meta_entries(sample_id, sample_meta_by_id):
            flat_data = dict(base_data)
            apply_mapping_fields(flat_data, meta)
            flat_data["max_cross_group_pairwise_iptm"] = compute_max_cross_group_pairwise_iptm(
                flat_data.get("pairwise_iptm"),
                flat_data.get("str_chain_map", "")
            )
            rows.append(flat_data)

    # Prune columns with no data across all rows
    writing_keys = prune_empty_columns(
        rows,
        writing_keys
    )

    # Write to CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=writing_keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(format_row_for_multiqc_csv(row, writing_keys))

def extract_alphafold3_results(json_files, samples, output_csv, sample_meta_by_id):
    rows = []
    all_keys = ["chain_iptm","chain_pair_iptm", "chain_pair_pae_min", "ptm", "iptm", "chain_ptm", "fraction_disordered", "has_clash", "ranking_score"]
    writing_keys = ["id", "model_input_id", "left_source_id", "right_source_id", "ptm", "iptm", "int_chain_map", "str_chain_map", "max_cross_group_chain_pair_iptm", "fraction_disordered", "has_clash", "ranking_score"]
    writing_keys_extra = set()
    # First pass: read and collect keys
    cntr = 0
    for file in json_files:
        with open(file) as f:
            data = json.load(f)
            sample_id = samples[cntr]
            cntr += 1
            for meta in get_meta_entries(sample_id, sample_meta_by_id):
                flat_data = {"id": sample_id}
                apply_mapping_fields(flat_data, meta)
                flat_data["max_cross_group_chain_pair_iptm"] = compute_max_cross_group_chain_pair_iptm(
                    data.get("chain_pair_iptm"),
                    flat_data.get("int_chain_map", "")
                )
                for key, value in data.items():
                    if key not in all_keys:
                        continue
                    if isinstance(value, dict) or isinstance(value, list):
                        sub_data = flatten_deep(value, index_lists=True)
                        for sub_key, sub_value in sub_data.items():
                            flat_data[f"{key}_{sub_key}"] = sub_value
                            writing_keys_extra.add(f"{key}_{sub_key}")
                    else:
                        flat_data[key] = value
                rows.append(flat_data)

    # Prune columns with no data across all rows
    final_writing_vals = writing_keys + sorted(writing_keys_extra)
    final_writing_vals = prune_empty_columns(
        rows,
        final_writing_vals,
    )

    # Write to CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=final_writing_vals)
        writer.writeheader()
        for row in rows:
            writer.writerow(format_row_for_multiqc_csv(row, final_writing_vals))



def main():
    parser = argparse.ArgumentParser(description="Collect confidence scorer generating by different models..")
    parser.add_argument("-i", "--input", help="Path to input directory containing all json files")
    parser.add_argument("-o", "--output", default="output_msa", help="output CSV file")
    parser.add_argument("--model", default="boltz", help="model generating the input files can be one of [boltz, alphafold3, and colabfold]")
    parser.add_argument("--sample-metadata", default=None, help="JSON file with list of sample metadata maps including id/report_id/int_chain_map/str_chain_map")

    args = parser.parse_args()
    json_files = [str(f) for f in Path(args.input).glob("*.json") if f.is_file()]
    sample_meta_by_id = {}
    if args.sample_metadata:
        with open(args.sample_metadata, "r") as f:
            sample_metadata = json.load(f)
        sample_meta_by_id = build_sample_meta_by_id(sample_metadata)

    if args.model.lower() == "boltz":
        extract_boltz_results(json_files, args.output, sample_meta_by_id)
    elif args.model.lower() == "colabfold":
        samples = [Path(p).name.split("_scores_")[0] for p in json_files]
        extract_colabfold_results(json_files, samples, args.output, sample_meta_by_id)
    elif args.model.lower() == "alphafold3":
        samples = [Path(p).name.split("_summary_confidences.json")[0] for p in json_files]
        extract_alphafold3_results(json_files, samples, args.output, sample_meta_by_id)

    print(f"CSV file created at: {args.output}")

if __name__ == "__main__":
    main()
