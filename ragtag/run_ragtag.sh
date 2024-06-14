#!/bin/bash

#SBATCH --job-name=test_run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --time=0-08:00:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/James_sepi/test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/James_sepi/test_output.txt

# Load necessary modules and activate the Python environment
module purge
module load Anaconda3/2023.09-0
source activate $DATA/python3_12_2

# Directory paths
ref_epi_dir="/data/biol-micro-genomics/kell7366/James_sepi/ref_epi/"
assembly_dir="/data/biol-micro-genomics/kell7366/James_sepi/Assemblies/"
fasta_dir="/data/biol-micro-genomics/kell7366/James_sepi/ragtag/fasta/"
temp_dir="/data/biol-micro-genomics/kell7366/James_sepi/ragtag/temp/"

# Create output and temporary directories if they don't exist
mkdir -p "$fasta_dir"
mkdir -p "$temp_dir"

# Loop through each assembly genome
for assembly_genome in "$assembly_dir"/*.fasta; do
  assembly_name=$(basename "$assembly_genome" | cut -d'_' -f1)
  final_scaffold_file="${fasta_dir}/${assembly_name}_scaffold.fasta"
  
  # Check if the final scaffold file already exists
  if [ -f "$final_scaffold_file" ]; then
    echo "Final scaffold file already exists for assembly $assembly_name. Skipping."
    continue
  fi

  input_genome="$assembly_genome"

  # Create a temporary directory for the current assembly
  temp_assembly_dir="${temp_dir}/${assembly_name}"
  mkdir -p "$temp_assembly_dir"

  # Loop through each reference genome in the ref_epi directory
  for ref_genome in "$ref_epi_dir"/*.fna; do
    ref_name=$(basename "$ref_genome" .fna)
    output_dir="${temp_assembly_dir}/${ref_name}"

    # Run ragtag scaffolding
    ragtag.py scaffold -o "$output_dir" "$ref_genome" "$input_genome"

    # Check if the ragtag scaffold file exists
    scaffold_file="${output_dir}/ragtag.scaffold.fasta"
    if [ -f "$scaffold_file" ]; then
      # Update input genome for the next iteration
      input_genome="$scaffold_file"
    else
      echo "File not found in directory $output_dir"
      break  # Exit the loop if the scaffold file is not found
    fi
  done

  # Move the final scaffold file to the fasta directory with the final name
  if [ -f "$input_genome" ]; then
    mv "$input_genome" "$final_scaffold_file"
  else
    echo "Final scaffold file not found for assembly $assembly_name"
  fi

  # Remove the temporary directory for the current assembly
  rm -rf "$temp_assembly_dir"
done

# Remove the temporary directory if it's empty
rmdir "$temp_dir" 2>/dev/null

