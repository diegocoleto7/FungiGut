#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.host_filtering    = true
params.bac_filtering     = true
params.out_dir           = "${launchDir}/preprocessing"
params.human_index       = "${launchDir}/resources/Bowtie2_Indexes/Host"
params.bac_index         = "${launchDir}/resources/Bowtie2_Indexes/UHGG"
params.cpus              = 6
params.single_end        = false
params.data_dir          = "data"

params.min_length        = 50
params.quality_threshold = 20

workflow {

    if (params.single_end) {
        Channel
            .fromPath(["${params.data_dir}/*.fastq.gz","${params.data_dir}/*.fq.gz"])
            .map { file -> tuple(file.baseName.replaceAll(/\.fastq(\.gz)?$/, ''), file) }
            .set{raw}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/*_{1,2}.fastq.gz",
                "${params.data_dir}/*_R{1,2}.fastq.gz",
                "${params.data_dir}/*_R{1,2}_001.fq.gz",
                "${params.data_dir}/*_{1,2}.fq.gz",
                "${params.data_dir}/*_R{1,2}.fq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
            .set{raw}
    }


    trimmed = qc_trimming(raw)

    host_cleaned = params.host_filtering ? host_filtering(trimmed.reads) : trimmed.reads
    final_reads  = params.bac_filtering  ? bac_filtering(host_cleaned.reads) : host_cleaned
}




process qc_trimming {
    cpus 6
    memory "6 GB"
    publishDir "${params.out_dir}/qc_trimmed", mode: 'link'
    tag "$id"
    
    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_filtered_R*.fastq.gz"), emit:reads
    tuple val(id), path("${id}_fastp.json"), emit:json
    tuple val(id), path("${id}.html"), emit:html


    script:
    if (params.single_end)
        """
        fastp -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --cut_front -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --cut_front -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """
}



process host_filtering {
    cpus params.cpus
    memory "8 GB"
    publishDir "${params.out_dir}/host_filtered", mode: 'link', subDir: true
    tag "$id"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_hostclean_R*.fastq.gz"), emit: reads
    tuple val(id), path("log/${id}_hostclean.log"), emit: log

    script:

    def log_out = "log/${id}_hostclean.log"

    if (params.single_end)
        """
        mkdir -p "log"

        bowtie2 -x ${params.human_index}/human \
            -U ${reads[0]} \
            --un-gz ${id}_hostclean_R1.fastq.gz \
            -p ${task.cpus} -S /dev/null --met-file $log_out
        """
    else
        """
        mkdir -p "log"

        bowtie2 -x ${params.human_index}/human \
            -1 ${reads[0]} -2 ${reads[1]} \
            --un-conc-gz ${id}_hostclean_R%.fastq.gz \
            -p ${task.cpus} -S /dev/null --met-file $log_out
        """
}


process bac_filtering {
    cpus params.cpus
    memory "8 GB"
    publishDir "${params.out_dir}/bac_filtered", mode: 'link', subDir: true
    tag "$id"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_bacclean_R*.fastq.gz"), emit: reads
    tuple val(id), path("log/${id}_bacclean.log"), emit: log


    script:
    
    def log_out = "log/${id}_bacclean.log"

    if (params.single_end)
        """
        mkdir -p "log"

        bowtie2 -x ${params.bac_index}/library \
            -U ${reads[0]} \
            --un ${id}_bacclean_R1.fastq.gz \
            -p ${task.cpus} -S /dev/null --met-file $log_out
        """
    else
        """
        mkdir -p "log"

        bowtie2 -x ${params.bac_index}/library \
            -1 ${reads[0]} -2 ${reads[1]} \
            --un-conc-gz ${id}_bacclean_R%.fastq.gz \
            -p ${task.cpus} -S /dev/null --met-file $log_out
        """
}
