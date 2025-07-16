#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
import json
import argparse
import re
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output_dir", required=True)
parser.add_argument("--report", required=True)
args = parser.parse_args()

INPUT_FILE = args.input
OUTPUT_DIR = args.output_dir
REPORT_FILE = args.report

os.makedirs(OUTPUT_DIR, exist_ok=True)
report = []

def genome_rank(entry):
    source_raw = entry.get("source_database", "")
    source_score = 0 if source_raw == "SOURCE_DATABASE_REFSEQ" else 1

    category = entry.get("assembly_info", {}).get("refseq_category", "")
    category_score = 0 if category == "reference genome" else 1

    level = entry.get("assembly_info", {}).get("assembly_level", "").lower()
    level_priority = {
        "complete genome": 0,
        "chromosome": 1,
        "scaffold": 2,
        "contig": 3
    }
    level_score = level_priority.get(level, 4)

    scaffold_count = entry.get("assembly_stats", {}).get("number_of_scaffolds", 999999)
    scaffold_score = scaffold_count if isinstance(scaffold_count, int) else 999999

    return (source_score, category_score, level_score, scaffold_score)

with open(INPUT_FILE, "r") as f:
    species_list = [line.strip() for line in f if line.strip()]

for species in species_list:
    species_dir = os.path.join(OUTPUT_DIR, species.replace(" ", "_"))
    os.makedirs(species_dir, exist_ok=True)

    try:
        summary_cmd = [
            "datasets", "summary", "genome", "taxon", species,
            "--as-json-lines", "--assembly-source", "all"
        ]
        result = subprocess.run(summary_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        lines = result.stdout.decode().strip().split("\n")

        genomes = [json.loads(line) for line in lines if line.strip()]
        if not genomes:
            report.append({"species": species, "accession": "NOT_FOUND"})
            continue

        best = sorted(genomes, key=genome_rank)[0]
        best_accession = best["accession"]

        zip_path = os.path.join(species_dir, "dataset.zip")
        download_cmd = [
            "datasets", "download", "genome", "accession", best_accession,
            "--filename", zip_path
        ]
        subprocess.run(download_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(["unzip", "-o", zip_path, "-d", species_dir], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        report.append({
            "species": species,
            "taxid": best.get("organism", {}).get("tax_id", "NA"),
            "accession": best_accession,
            "source_database": best.get("source_database", ""),
            "assembly_level": best.get("assembly_info", {}).get("assembly_level", ""),
            "refseq_category": best.get("assembly_info", {}).get("refseq_category", ""),
            "scaffold_count": best.get("assembly_stats", {}).get("number_of_scaffolds", "NA")
        })

    except subprocess.CalledProcessError as e:
        report.append({
            "species": species,
            "taxid": "ERROR",
            "accession": "ERROR",
            "source_database": "ERROR",
            "assembly_level": "ERROR",
            "refseq_category": "ERROR",
            "scaffold_count": "ERROR",
            "error": str(e)
        })

for dir_name in os.listdir(OUTPUT_DIR):
    old_path = os.path.join(OUTPUT_DIR, dir_name)

    if not os.path.isdir(old_path):
        continue

    clean_name = re.sub(r"[\[\]]", "", dir_name)
    clean_path = os.path.join(OUTPUT_DIR, clean_name)

    if clean_name != dir_name:
        if not os.path.exists(clean_path):
            os.rename(old_path, clean_path)
        else:
            shutil.rmtree(old_path)  
        old_path = clean_path

    for root, _, files in os.walk(old_path):
        for file_name in files:
            if file_name.endswith(".fna"):
                src_fna = os.path.join(root, file_name)
                dst_name = f"{clean_name}_{file_name}"
                dst_fna = os.path.join(OUTPUT_DIR, dst_name)
                shutil.move(src_fna, dst_fna)

    shutil.rmtree(old_path)


df = pd.DataFrame(report)
df.to_csv(REPORT_FILE, index=False)
