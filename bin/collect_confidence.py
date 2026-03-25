#!/usr/bin/env python3

import os, sys
import json
import csv
import statistics
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

def extract_alphafold3_results(json_files, samples, output_csv):
    rows = []
    all_keys = ["chain_iptm","chain_pair_iptm", "chain_pair_pae_min", "ptm", "iptm", "chain_ptm", "fraction_disordered", "has_clash", "ranking_score"]
    writing_keys = ["id", "ptm", "iptm", "fraction_disordered", "has_clash", "ranking_score"]
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
                if isinstance(value, dict) or isinstance(value, list):
                    sub_data = flatten_deep(value, index_lists=True)
                    for sub_key, sub_value in sub_data.items():
                        flat_data[f"{key}_{sub_key}"] = sub_value
                        writing_keys_extra.add(f"{key}_{sub_key}")
                else:
                    flat_data[key] = value
            rows.append(flat_data)

    # Write to CSV
    final_writing_vals = writing_keys + sorted(writing_keys_extra)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=final_writing_vals)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in final_writing_vals})



def main():
    parser = argparse.ArgumentParser(description="Collect confidence scorer generating by different models..")
    parser.add_argument("-i", "--input", help="Path to input directory containing all json files")
    parser.add_argument("-o", "--output", default="output_msa", help="output CSV file")
    parser.add_argument("--model", default="boltz", help="model generating the input files can be one of [boltz, alphafold3, and colabfold]")

    args = parser.parse_args()
    json_files = [str(f) for f in Path(args.input).glob("*.json") if f.is_file()]

    if args.model.lower() == "boltz":
        samples = [Path(p).name.split("_boltz_")[0].split("confidence_")[1] for p in json_files]
        extract_boltz_results(json_files, samples, args.output)
    elif args.model.lower() == "colabfold":
        samples = [Path(p).name.split("_scores_")[0] for p in json_files]
        extract_colabfold_results(json_files, samples, args.output)
    elif args.model.lower() == "alphafold3":
        samples = [Path(p).name.split("_summary_confidences.json")[0] for p in json_files]
        extract_alphafold3_results(json_files, samples, args.output)

    print(f"CSV file created at: {args.output}")

if __name__ == "__main__":
    main()
