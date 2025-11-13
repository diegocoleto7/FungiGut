nextflow.enable.dsl=2

params.out_dir              = "resources"
params.genome_list          = "assets/species_list.txt"
params.download_host        = true
params.download_prokaryote  = true
params.update_db            = false
params.make_accession_info  = true
params.threads              = 8

workflow {

    if (params.download_host) {
        host_index = download_human_index()
    }

    if (params.download_prokaryote) {
        bact_index = download_prokaryote_index()
    }

    if (params.update_db) {
        fungi = download_fungi_genomes(file(params.genome_list))

        index_fungi(fungi.genomes)
            .filtered_fasta
            .set { fungi_fna }

        if (params.make_accession_info) {
            make_accession_info(fungi_fna)
        }

    } else {
        fungigutdb = download_fungigutdb()

    }
}

process download_human_index {
    publishDir "${params.out_dir}/Bowtie2_Indexes", mode: 'link'

    output:
    path("Host")

    script:
    """
    aria2c -x 8 -s 8 -o Host.tar.gz https://zenodo.org/records/17581472/files/Host.tar.gz?download=1
    tar -xzf Host.tar.gz
    rm -f Host.tar.gz
    """
}


process download_prokaryote_index {
    publishDir "${params.out_dir}/Bowtie2_Indexes", mode: 'link'

    output:
    path("UHGG")

    script:
    """
    aria2c -x 8 -s 8 -o UHGG.tar.gz https://zenodo.org/records/17581472/files/UHGG.tar.gz?download=1
    tar -xzf UHGG.tar.gz
    rm -f UHGG.tar.gz    
    """
}

process download_fungigutdb {
    publishDir "${params.out_dir}", mode: 'link'

    output:
    path("FungiGut_db")

    script:
    """
    aria2c -x 8 -s 8 -o FungiGut_db.tar.gz https://zenodo.org/records/17581472/files/FungiGut_db.tar.gz?download=1
    tar -xzf FungiGut_db.tar.gz
    rm -f FungiGut_db.tar.gz
    """
}

process download_fungi_genomes {
    publishDir "${params.out_dir}/FungiGut_db", mode: 'link'

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
    publishDir "${params.out_dir}/FungiGut_db", mode: 'link'

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

    publishDir "${params.out_dir}/FungiGut_db", mode: 'link'

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
