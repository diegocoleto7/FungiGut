# FungiGut 🍄

This repository provides a suite of **Nextflow-based workflows** to build genome databases, preprocess samples, and profile fungal abundance from shotgun metagenomic (WGS) reads of the human gut microbiome.

---

## 📂 Repository Structure

```
FungiGut/
├── resources.nf        # Workflow for downloading and indexing host, bacteria and fungal genomes.
├── preprocess.nf       # Workflow for quality control filtering and host/bac reads removal.
├── taxprofiler.nf      # Workflow for alignment and fungal abundance profiling.
├── bin/                # Auxiliary scripts
├── assets/             # Fungi species list
├── FungiGut.yml        # Conda environment
└── README.md   
```

- **resources.nf**  
  - Downloads genomes (human, bacteria, fungi).
  - Filter quality on fungi genomes.  
  - Builds indices for Bowtie2 and BWA.  

- **preprocess.nf**  
  - Runs Fastp for read quality assessment and trimming.  
  - Filters out host and bacterial reads using prebuilt indices.  

- **taxprofiler.nf**  
  - Aligns cleaned reads against the fungal reference database.  
  - Calculates relative abundances (e.g., TPM, RPKM) with tools like CoverM or MetaPhlAn.  

- **bin/**  
  - Helper scripts for results, generate summary reports, etc.  

- **assets/**  
  - Contains FungiGut species list.  

- **FungiGut.yml**  
  - Defines a Conda environment.

---

## ⚙️ Installation & Setup

### 1. Prerequisites

- **Java 8+** (required by Nextflow)  
- **Conda** or **Miniconda** (for dependency management)  
- **Git** (to clone the repository)

### 2. Clone the Repository

```bash
git clone https://github.com/diegocoleto7/FungiGut.git
cd FungiGut
```
### 3. Create & Activate the Conda Environment

```bash
conda env create -n fungigut -f FungiGut.yml
conda activate fungigut
```
### 4. (Optional) Define a Central Project Directory

To keep multiple projects organized, you can set an environment variable pointing to the FungiGut scripts:
```bash
export FUNGI=/path/to/FungiGut
```
Then in any project folder:
```bash
cd /my/project/folder
nextflow run $FUNGI/resources.nf   # or preprocess.nf, taxprofiler.nf...
```


| Workflow         | Parameter             | Default                                | Description & Tips                                                  |
| ---------------- | --------------------- | -------------------------------------- | ------------------------------------------------------------------- |
| **resources.nf** | `out_dir`             | `${launchDir}/resources`               | Output directory. Change if you need a different location.          |
|                  | `genome_list`         | `/assets/species_list.txt`             | Path to fungal species list. Edit to add or remove species.         |
|                  | `download_host`       | `true`                                 | Download human genome. Set to `false` if already available.         |
|                  | `download_bacteria`   | `true`                                 | Download bacterial database (UHGG).                                 |
|                  | `download_fungi`      | `true`                                 | Download fungal genomes.                                            |
|                  | `make_accession_info` | `false`                                | Generate accession→taxid mapping. Increases runtime and disk usage. |
|                  | `maxcpus`             | `6`                                    | Maximum threads for indexing. Adjust to your hardware.              |
