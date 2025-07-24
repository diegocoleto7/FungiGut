nextflow.enable.dsl=2

params.out_dir              = "resources"
params.genome_list          = "assets/species_list.txt"
params.download_host        = true
params.download_bacteria    = true
params.download_fungi       = true
params.make_accession_info  = false
params.threads           = 10

workflow {

    if (params.download_host) {
        host_fasta = download_human_genome()
        index_human_bowtie2(host_fasta)
    }

    if (params.download_bacteria) {
        bact_fasta = download_bacteria_genome()
        index_bact_bowtie2(bact_fasta)
    }

    if (params.download_fungi) {
        fungi = download_fungi_genomes(file(params.genome_list))

        index_fungi(fungi.genomes)
            .filtered_fasta
            .set { fungi_fna }

        if (params.make_accession_info) {
            make_accession_info(fungi_fna)
        }
    }
}

process download_human_genome {
    cpus 2
    publishDir "${params.out_dir}/Host", mode: 'copy'

    output:
    path("human.fna")

    script:
    """
    aria2c -x 8 -s 8 -o human.fna.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
    pigz -d -f human.fna.gz

    """
}

process index_human_bowtie2 {
    cpus params.threads
    publishDir "${params.out_dir}/Bowtie2_Indexes/Host", mode: 'copy'

    input:
    path fasta_file

    output:
    path("${fasta_file.simpleName}.*")

    script:
    """
    bowtie2-build --threads ${task.cpus} "${fasta_file}" "${fasta_file.simpleName}"
    """
}

process download_bacteria_genome {
    publishDir "${params.out_dir}/UHGG", mode: 'copy'

    output:
    path("library.fna")

    script:
    """
    aria2c -x 8 -s 8 -o library.fna https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/kraken2_db_uhgg_v2.0.2/library/library.fna

    """
}

process index_bact_bowtie2 {
    cpus params.threads
    publishDir "${params.out_dir}/Bowtie2_Indexes/UHGG", mode: 'copy'

    input:
    path fasta_file

    output:
    path("${fasta_file.simpleName}.*")

    script:
    """
    bowtie2-build --threads ${task.cpus} "${fasta_file}" "${fasta_file.simpleName}"
    """
}

process download_fungi_genomes {
    publishDir "${params.out_dir}/FungiGut_db", mode: 'copy'

    input:
    path species_list

    output:
    path("genomes"), emit: genomes
    path("genome_report.csv"), emit: report

    script:
    """
    mkdir -p genomes
    download_fungi_genomes.py \
        --input ${species_list} \
        --output_dir genomes \
        --report genome_report.csv
    """
}

process index_fungi {
    publishDir "${params.out_dir}/FungiGut_db", mode: 'copy'

    input:
    path genome_dir

    output:
    path("genomes.fna"), emit: filtered_fasta
    path("FungiGut.*")

    script:
    """
    cat ${genome_dir}/*.fna | seqkit seq -m 500 > genomes.fna
    bwa index -a bwtsw -p FungiGut genomes.fna
    """
}

process make_accession_info {
    cpus 4
    memory '4 GB'

    publishDir "${params.out_dir}/FungiGut_db", mode: 'copy'

    input:
    path(fasta_file)

    output:
    path("accession2info.txt")

    script:
    """
    get_len.py ${fasta_file} > length.tsv

    aria2c -x 4 -s 4 -o nucl_gb.accession2taxid.gz \
        https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    aria2c -x 4 -s 4 -o nucl_wgs.accession2taxid.gz \
        https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz

    pigz -d -f nucl_gb.accession2taxid.gz
    pigz -d -f nucl_wgs.accession2taxid.gz

    cat nucl_gb.accession2taxid nucl_wgs.accession2taxid > all_accessions.tsv
    awk -F'\t' '{ print \$2 "\t" \$3 }' all_accessions.tsv > acc__taxid.tsv

    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 \
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

    mkdir taxdump
    tar -xzf taxdump.tar.gz --directory taxdump names.dmp nodes.dmp

    generate_accession2info.py length.tsv acc__taxid.tsv taxdump/nodes.dmp taxdump/names.dmp accession2info.txt
    """
}
