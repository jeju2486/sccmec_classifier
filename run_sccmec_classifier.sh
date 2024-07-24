#!/bin/bash

# Set the path to your query sequence
query_sequence="/data/biol-micro-genomics/kell7366/mGWAS_staphy/sccmec_classifier/gene_db.fasta"

# Set the directory containing reference sequences and flanking length
reference_dir="/data/biol-micro-genomics/kell7366/sccmec_ncbi_dataset/filtered_data/ref_fas/ref_aureus/all"

# Create a directory to store the files
output_dir="/data/biol-micro-genomics/kell7366/mGWAS_staphy/sccmec_classifier"

save_temp="TRUE"

module purge
module load Anaconda3/2023.09-0
source activate $DATA/sccmec_classifier_env

bash run_script.sh -q "$query_sequence" -r "$reference_dir" -o "$output_dir" -s "$save_temp"

python merge.py