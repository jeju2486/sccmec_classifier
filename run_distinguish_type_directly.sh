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
query_sequence="/data/biol-micro-genomics/kell7366/James_sepi/mec_class.fasta"

# Set the directory containing reference sequences and flanking length
reference_dir="/data/biol-micro-genomics/kell7366/James_sepi/ragtag/fasta"

# Create a directory to store the files
output_dir="/data/biol-micro-genomics/kell7366/James_sepi"

mkdir -p "$output_dir/bed_direct"
mkdir -p "$output_dir/sam_direct"
mkdir -p "$output_dir/filtered_sam_direct"

#Iterate through files in the reference directory
for reference_genome in "$reference_dir"/*.fasta; do
    # Get the reference sequence name (without extension) to use in the output SAM file name
    reference_name=$(basename "$reference_genome" | cut -d '.' -f 1)
    if [ -e "$output_dir/sam_direct/${reference_name}.sam" ]; then
      echo "File $file already exists. Skipping."
    
    else
        # Your processing logic here
        echo "Processing $reference_name"

    # Run Minimap2 and save the SAM file in the SAM directory
    minimap2 -a "$reference_genome" "$query_sequence" > "$output_dir/sam_direct/${reference_name}.sam"
    fi
done

#loda python module
module purge
module load Anaconda3/2023.09-0
conda activate $DATA/python3_12_2

#filter sam file high match rate with enough lengh
python sam_filter.py

#load modules
module purge
module load BEDTools/2.30.0-GCC-11.2.0


for reference_genome in "$reference_dir"/*.fasta; do
    # Get the reference sequence name (without extension) to use in the output SAM file name
    reference_name=$(basename "$reference_genome" | cut -d '.' -f 1)
    
    echo "Processing $reference_name"
    
    # Convert BAM file into a BED file and save it in the BED directory
    bedtools bamtobed -i "$output_dir/filtered_sam_direct/${reference_name}.sam" > "$output_dir/bed_direct/${reference_name}.bed"

    # Create genome length information file and save it in the BED directory
    # The score data is calculated by MAPQ: MAPping Quality. It equals -10 log10 Pr{mapping position is wrong}
    cat "$reference_genome" | awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' | sed 's/>//; s/ / /' > "$output_dir/bed_direct/${reference_name}_genome_size.txt"
    
done

#loda python module
module purge
module load Anaconda3/2023.09-0
conda activate $DATA/python3_12_2

#extracting mec type according to type info
python extract_mectype.py

#extracting ccr type by best match
python extract_ccr.py

#calculate the sccmec type corresponding ccr and mec type
python combining_files.py