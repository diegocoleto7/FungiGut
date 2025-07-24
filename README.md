# FungiGut

| Workflow         | Parameter             | Default                                | Description & Tips                                                  |
| ---------------- | --------------------- | -------------------------------------- | ------------------------------------------------------------------- |
| **resources.nf** | `out_dir`             | `${launchDir}/resources`               | Output directory. Change if you need a different location.          |
|                  | `genome_list`         | `${scriptDir}/assets/species_list.txt` | Path to fungal species list. Edit to add or remove species.         |
|                  | `download_host`       | `true`                                 | Download human genome. Set to `false` if already available.         |
|                  | `download_bacteria`   | `true`                                 | Download bacterial database (UHGG).                                 |
|                  | `download_fungi`      | `true`                                 | Download fungal genomes.                                            |
|                  | `make_accession_info` | `false`                                | Generate accession→taxid mapping. Increases runtime and disk usage. |
|                  | `maxcpus`             | `6`                                    | Maximum threads for indexing. Adjust to your hardware.              |
