#!/bin/bash

# Set the path to your query sequence
query_sequence="/data/biol-micro-genomics/kell7366/sccmec_classifier/gene_db.fasta"

# Set the directory containing reference sequences and flanking length
reference_dir="/data/biol-micro-genomics/kell7366/sccmec_classifier/Assemblies/"

# Create a directory to store the files
output_dir="/data/biol-micro-genomics/kell7366/sccmec_classifier"

save_temp="TRUE"

module purge
module load Anaconda3/2023.09-0

python setup.py

bash run_script.sh -q "$query_sequence" -r "$reference_dir" -o "$output_dir" -s "$save_temp"
