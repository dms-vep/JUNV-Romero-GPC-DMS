# Description:
# Downloads accessions and process/extract sequences and metadata

# Author:
# Caleb Carr

rule download_and_process_accessions:
    """
    This rule runs the download_NCBI_sequences.py script which
    downloads all genbank files from a list of accessions and 
    extracts metadata. 
    """
    input:
        accession_list = config["Accession_list"],
    params:
        local_file_dir = config["Local_genbanks"],
        genome_size_threshold_lower = lambda wildcards: config["Genome_size_threshold_lower"][wildcards.gene],
        genome_size_threshold_upper = lambda wildcards: config["Genome_size_threshold_upper"][wildcards.gene],
        max_frac_N = lambda wildcards: config["max_frac_N"][wildcards.gene],
        accesstions_to_exclude = lambda wildcards: config["Accessions_to_exclude"][wildcards.gene],
        accessions_to_include = lambda wildcards: config["Accessions_to_include"][wildcards.gene],
        remove_duplicates = config["Remove_duplicates"],
        asia = config["Asia"],
        oceania = config["Oceania"],
        africa = config["Africa"],
        europe = config["Europe"],
        south_america = config["South America"],
        north_america = config["North America"],
        genbank_features = config["Genbank_features"],
        host_grouper = config["Host_grouper"],
        virus_grouper = lambda wildcards: config["Virus_groupings"][wildcards.gene],
        gene_names = lambda wildcards: config["Gene_groupings"][wildcards.gene],
        desired_gene = lambda wildcards: wildcards.gene,
        reference_genome = lambda wildcards: config["Reference_fasta"][wildcards.gene],
        temp_fasta_file = config["Temp_fasta_file"],
        temp_alignment_file = config["Temp_alignment_file"],
        min_id_cutoff = lambda wildcards: config["Min_id_cutoff"][wildcards.gene],
    output:
        fasta_sequences = "Results/{gene}/nucleotide.fasta",
        metadata = "Results/{gene}/metadata.tsv",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/download_and_process_accessions.txt",
    script:
        "../Scripts/download_NCBI_sequences.py"


rule create_new_references:
    """
    This rule runs a script to create new references.
    """
    input:
        lambda wildcards: config["Reference_genbank"][wildcards.gene],
    params:
        gene_name = lambda wildcards: config["Gene_groupings"][wildcards.gene][0],
    output:
        new_reference = "Results/{gene}/{gene}_reference.gb",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/create_new_references.txt",
    shell:
        "python Scripts/newreference.py "
        "--reference {input} "
        "--output-genbank {output.new_reference} "
        "--new-name GPC "
        "--gene {params.gene_name} &> {log}"
