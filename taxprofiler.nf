#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.single_end         = false
params.reads_dir          = "${launchDir}/preprocessing/bac_filtered"
params.fungi_dir          = "${launchDir}/resources/FungiGut_db"
params.results_dir        = "${launchDir}/results"
params.accession_info     = "${launchDir}/resources/FungiGut_db/accession2info.txt"
params.out_rds            = "phyloseq_object.rds"
params.min_map            = 100
params.max_ed             = 5
params.pct_id             = 0.97
params.read_cutoff        = 1
params.cpus               = 6

params.run_alignment      = true

workflow {
    reads = params.single_end ?
        Channel.fromPath("${params.reads_dir}/*_R1.fastq.gz")
               .map { f -> tuple(f.baseName.replaceAll('_R1',''), [f]) }
        :
        Channel.fromFilePairs([
                        "${params.reads_dir}/*_R{1,2}.fastq.gz",
                        "${params.reads_dir}/*_R{1,2}.fq.gz",
                ])
               .ifEmpty { error "Cannot find reads in ${params.reads_dir}!" }

    sam_files = params.run_alignment ?
        bwa_align(reads)
        :
        Channel
          .fromPath("${params.results_dir}/sam/*.sam")
          .map { f -> tuple(f.simpleName, f) }


    abundance = compute_abundance(sam_files)

    done_signal = abundance
        .map { id, txt -> txt }
        .collect()
        .map { "${params.results_dir}/abundance_txt" }

    rds = generate_rds(done_signal)
}

process bwa_align {
    cpus params.cpus
    publishDir "${params.results_dir}/sam", mode: 'copy'
    tag "$id"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.sam"), emit:sam

    script:
    if (params.single_end)
        """
        bwa mem -t ${task.cpus} -a "${params.fungi_dir}/FungiGut" ${reads[0]} > ${id}.sam
        """
    else
        """
        bwa mem -t ${task.cpus} -a "${params.fungi_dir}/FungiGut" ${reads[0]} ${reads[1]} > ${id}.sam
        """
}

process compute_abundance {
    cpus 2
    memory "4 GB"
    publishDir "${params.results_dir}/abundance_txt", mode: 'copy'
    tag "$id"

    input:
    tuple val(id), path(sam_file)

    output:
    tuple val(id), path("${id}.txt")

    script:
    """
    compute-abundances.py \\
        ${sam_file} \\
        --accession_info ${params.accession_info} \\
        --output ${id}.txt \\
        --min_map ${params.min_map} \\
        --max_ed ${params.max_ed} \\
        --pct_id ${params.pct_id} \\
        --read_cutoff ${params.read_cutoff} \\
        --raw_counts
    """
}

process generate_rds {
    cpus 2
    publishDir "${params.results_dir}/phyloseq", mode: 'copy'

    input:
    val abundance_dir

    output:
    path("*.rds")

    script:
    """
    abundancetophyseq.R ${abundance_dir} ${params.out_rds}
    """
}
