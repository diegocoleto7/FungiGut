# FungiGut 

This repository provides a suite of **Nextflow-based workflows** to **build genome databases**, **preprocess samples**, and **profile fungal abundance** from **shotgun metagenomic** (WGS) reads of the human gut microbiome.

---

##  Repository Structure

```
FungiGut/
├── resources.nf        # Download prebuilt Host/UHGG Bowtie2 indices; download/build FungiGutDB (BWA + accession2info)
├── preprocess.nf       # fastp QC+trimming; host and bacterial filtering with Bowtie2
├── taxprofiler.nf      # BWA-MEM alignment; MiCoP-like abundances; Phyloseq .rds
├── bin/                # Auxiliary scripts (e.g., compute-abundances.py, R utilities)
├── assets/             # species_list.txt (fungal species list)
├── FungiGut.yml        # Conda environment
└── README.md
```
---

##  Installation & Setup

#### Prerequisites

- **[Java 17 (or later,up to 24)](http://www.oracle.com/technetwork/java/javase/downloads/index.html)** (required by Nextflow)
- **Conda** or **Miniconda** (for dependency management) 
- **Git** (to clone the repository)

#### 1) Clone the Repository

```bash
git clone https://github.com/diegocoleto7/FungiGut.git
cd FungiGut
```


#### 2) Create & Activate the Conda Environment

```bash
conda env create -n fungigut -f FungiGut.yml
conda activate fungigut
```



#### 3) (Optional) Define a Central Project Directory

To keep multiple projects organized, you can set an environment variable pointing to the FungiGut scripts:
```bash
export FUNGI=/path/to/FungiGut
```
Then in any project folder:
```bash
cd /my/project/folder
nextflow run $FUNGI/resources.nf   # or preprocess.nf, taxprofiler.nf...
```




##  Usage Overview


| Workflow           | Output Directory | Description                              |
| ------------------ | ---------------- | ---------------------------------------- |
| **resources.nf**   | `resources/`     | Host/UHGG indices + FungiGutDB           |
| **preprocess.nf**  | `preprocessing/` | QC, trimming, host & bacterial filtering |
| **taxprofiler.nf** | `results/`       | BWA alignments, abundances, Phyloseq RDS |



---
### 1) Resources — `resources.nf`



 **Default Parameters**
```groovy
params.out_dir             = "resources"
params.genome_list         = "assets/species_list.txt"
params.download_host       = true
params.download_prokaryote = true
params.update_db           = false
params.make_accession_info = true
params.threads             = 8

```




 **What It Does**

* **Download prebuilt Bowtie2 indices:**
    * **Host** → `--download_host true`.
    * **UHGG** → `--download_prokaryote true`.
* **Either** build **FungiGutDB** from assets/species_list.txt`--update_db true`,
* **Or** download FungiGut_db from Zenodo **(recomended)**--update_db false.




>**Note**:`make_accession_info` runs only when you rebuild (it relies on the FASTA you just generated).



---
 **Example Runs**
##### **A) Default: download Host + UHGG + FungiGut_db from Zenodo**
```bash
nextflow run $FUNGI/resources.nf
```
##### **B) Rebuild FungiGutDB from a species list**
```bash
nextflow run $FUNGI/resources.nf \
  --update_db true
```

### 2) Preprocessing — `preprocess.nf`

**Default parameters**
```groovy
params.host_filtering    = true
params.bac_filtering     = true
params.out_dir           = "${launchDir}/preprocessing"
params.human_index       = "${launchDir}/resources/Bowtie2_Indexes/Host"
params.bac_index         = "${launchDir}/resources/Bowtie2_Indexes/UHGG"
params.cpus              = 8
params.single_end        = false
params.data_dir          = "data"

params.trim_front        = 5
params.min_length        = 50
params.quality_threshold = 20
```



**What It Does**
* Performs **QC and trimming** with *fastp* (JSON + HTML reports).
* Applies **host filtering** if `--host_filtering true`.
* Applies **bacterial filtering** if `--bac_filtering true`.
* Produces **clean** FASTQs:
    * After host filtering:`<ID>_hostclean_R1/R2.fastq.gz`
    * After bacteria filtering:`<ID>_bacclean_R1/R2.fastq.gz`
 
 
 
 
**Input File Patterns**

* **Paired-end**: any of `*_R{1,2}_001.fastq.gz`, `_{1,2}.fastq.gz`, `_R{1,2}.fastq.gz` (and `.fq.gz` variants).
* **Single-end**: `*.fastq.gz` or `*.fq.gz` with `--single_end true`.

---
**Example Runs**

##### **A) Default (Paired-end with host + bacteria filtering)**
```bash
nextflow run $FUNGI/preprocess.nf \
  --data_dir /path/to/raw/reads \
  --cpus 12
```

##### **B) Single-end**
```bash
nextflow run $FUNGI/preprocess.nf \
  --data_dir /path/to/single/end/raw/reads \
  --single_end true
```

### 3) Fungal Profiling — `taxprofiler.nf`
**Default Parameters**
```groovy
params.single_end     = false
params.reads_dir      = "${launchDir}/preprocessing/bac_filtered"
params.fungi_dir      = "${launchDir}/resources/FungiGut_db"
params.results_dir    = "${launchDir}/results"
params.accession_info = "${launchDir}/resources/FungiGut_db/accession2info.txt"
params.out_rds        = "phyloseq_object.rds"

# compute-abundances filters
params.min_map        = 100
params.max_ed         = 1
params.pct_id         = 0.99
params.read_cutoff    = 10

params.cpus           = 6
params.raw_counts     = false
params.run_alignment  = true

```
**What It Does**

* Performs **BWA-MEM alignment** against FungiGutDB `--run_alignment true`.
* Estimates **abundance** (MiCoP) with filters:
    * `min_map` , `max_ed`, `pct_id`, `read_cutoff`.
* Generates a **Phyloseq object** (`.rds`) from sample abundance tables.

**Expected Inputs**
* By default, FASTQs from `preprocessing/bac_filtered`:
    * Paired-end: `*_R{1,2}.fastq.gz`.
    * Single-end: `*_R1.fastq.gz` with `--single_end true`.
* If skipping alignment (`--run_alignment false`): place SAMs in `results/sam/` and set `--results_dir` accordingly.
---

**Example Runs**

##### **A) Default: Full run (align + abundances + RDS)**
```bash
nextflow run $FUNGI/taxprofiler.nf \
  --cpus 12 \
```
##### **B) Use existing SAMs (skip alignment)**

```bash
nextflow run $FUNGI/taxprofiler.nf \
  --run_alignment false \
  --cpus 12
```
##### **C) Adjust assignment filters**
```bash
nextflow run $FUNGI/taxprofiler.nf \
  --min_map <> --max_ed <> --pct_id <> --read_cutoff <> \
  --cpus 12
```
