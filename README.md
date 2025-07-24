# FungiGut
This repository contains a set of Nextflow-based workflows for building genome databases, preprocessing samples and profiling WGS reads from the gut fungal microbiome.
```
FungiGut/
├── resources.nf        # Workflow for downloading and indexing host, bacteria and fungal genomes.
├── preprocess.nf       # Workflow for quality control filtering and host/bac reads removal.
├── taxprofiler.nf      # Workflow for alignment and fungal abundance profiling.
├── bin/                # Auxiliary scripts
├── assets/             # Example data and lists
├── FungiGut.yml        # Conda environment
└── README.md   
```
## Instalation
**1.- Clone the repository:**
```
git clone https://github.com/diegocoleto7/FungiGut.git
```
| Workflow         | Parameter             | Default                                | Description & Tips                                                  |
| ---------------- | --------------------- | -------------------------------------- | ------------------------------------------------------------------- |
| **resources.nf** | `out_dir`             | `${launchDir}/resources`               | Output directory. Change if you need a different location.          |
|                  | `genome_list`         | `${scriptDir}/assets/species_list.txt` | Path to fungal species list. Edit to add or remove species.         |
|                  | `download_host`       | `true`                                 | Download human genome. Set to `false` if already available.         |
|                  | `download_bacteria`   | `true`                                 | Download bacterial database (UHGG).                                 |
|                  | `download_fungi`      | `true`                                 | Download fungal genomes.                                            |
|                  | `make_accession_info` | `false`                                | Generate accession→taxid mapping. Increases runtime and disk usage. |
|                  | `maxcpus`             | `6`                                    | Maximum threads for indexing. Adjust to your hardware.              |
