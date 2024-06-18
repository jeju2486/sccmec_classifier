#!/bin/bash

#SBATCH --job-name=test_run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=devel
#SBATCH --time=0-00:03:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/James_sepi/test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/James_sepi/test_output.txt

#load modules
module purge
module load minimap2/2.24-GCCcore-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0

# Set the path to your query sequence
query_sequence="/data/biol-micro-genomics/kell7366/James_sepi/cassette_genes.fasta"

# Set the directory containing reference sequences and flanking length
reference_dir="/data/biol-micro-genomics/kell7366/James_sepi/ragtag/fasta"

# Create a directory to store the files
output_dir="/data/biol-micro-genomics/kell7366/James_sepi"

mkdir -p "$output_dir/bed"
mkdir -p "$output_dir/sam"
mkdir -p "$output_dir/temp"

#Iterate through files in the reference directory
for reference_genome in "$reference_dir"/*.fasta; do
    # Get the reference sequence name (without extension) to use in the output SAM file name
    reference_name=$(basename "$reference_genome" | cut -d '.' -f 1)
    if [ -e "$output_dir/sam/${reference_name}.sam" ]; then
      echo "File $file already exists. Skipping."
    
    else
        # Your processing logic here
        echo "Processing $reference_name"

    # Run Minimap2 and save the SAM file in the SAM directory
    minimap2 -x map-ont --secondary=no -a "$reference_genome" "$query_sequence" > "$output_dir/sam/${reference_name}.sam"

    # Convert BAM file into a BED file and save it in the BED directory
    bedtools bamtobed -i "$output_dir/sam/${reference_name}.sam" > "$output_dir/bed/${reference_name}.bed"

    # Create genome length information file and save it in the BED directory
    # The score data is calculated by MAPQ: MAPping Quality. It equals -10 log10 Pr{mapping position is wrong}
    cat "$reference_genome" | awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' | sed 's/>//; s/ / /' > "$output_dir/bed/${reference_name}_genome_size.txt"
    
    fi
done

#python code load
module purge
module load Anaconda3/2023.09-0
source activate $DATA/python3_12_2

python distinguish_mec_type.py
python extract_ccr.py

#calculate the sccmec type corresponding ccr and mec type
python combining_files.py
