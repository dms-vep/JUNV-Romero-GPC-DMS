#!/bin/bash

python gga_codon_muts_oligo_design.py \
    --tiles_csv Data/JUNV_GAA_assembly_fragments.csv \
    --mutations_to_make_csv Data/JUNV_mutations_to_make.csv \
    --avoid_motifs "TCTAGA" "CGTCTC" \
    --output_oligos_fasta Results/output_oligos.fa
