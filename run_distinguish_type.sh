#!/bin/bash

# Function to display usage message
usage() {
  echo "Usage: $0 -q <query_sequence> -r <reference_dir> -o <output_dir>"
  exit 1
}

# Parse command-line arguments
while getopts q:r:o: flag
do
    case "${flag}" in
        q) query_sequence=${OPTARG};;
        r) reference_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        \?) usage;;
    esac
done

mkdir -p "$output_dir/bed" "$output_dir/sam" "$output_dir/temp"

# Function to filter best matching result from SAM file
filter_best_match() {
  local sam_file=$1
  local filtered_sam_file=$2

  # Retain the SAM header and filter the alignments
  awk 'BEGIN {FS="\t"; OFS="\t"} /^@/ {print $0; next} !/^@/ {if (!best[$1] || $5 > best[$1]) {best[$1] = $5; align[$1] = $0}} END {for (read in align) print align[read]}' "$sam_file" > "$filtered_sam_file"
}

# Function to find and realign soft-clipped regions
extract_and_realign_soft_clipped() {
  local sam_file=$1
  local reference_genome=$2
  local output_dir=$3
  local reference_name=$4
  local output_suffix=$5

  # Extract soft-clipped sequences and realign
  awk 'BEGIN {FS="\t"; OFS="\t"} !/^@/ && $6 ~ /S/ {
    match($6, /[0-9]+S/)
    if (RSTART == 1) {
      clipped_length = substr($6, 1, RLENGTH - 1)
      clipped_seq = substr($10, 1, clipped_length)
    } else {
      clipped_length = substr($6, RSTART, RLENGTH - 1)
      clipped_seq = substr($10, length($10) - clipped_length + 1, clipped_length)
    }
    print ">"$1"\n"clipped_seq
  }' "$sam_file" > "$output_dir/temp/${reference_name}_soft_clipped${output_suffix}.fasta"

  minimap2 -a "$reference_genome" "$output_dir/temp/${reference_name}_soft_clipped${output_suffix}.fasta" > "$output_dir/temp/${reference_name}_soft_clipped${output_suffix}.sam"

  # Filter best match from the realigned soft-clipped SAM file
  filter_best_match "$output_dir/temp/${reference_name}_soft_clipped${output_suffix}.sam" "$output_dir/temp/${reference_name}_filtered_soft_clipped${output_suffix}.sam"
}

# Function to process the alignment and soft clipping
process_alignment() {
  local query_sequence=$1
  local reference_genome=$2
  local output_dir=$3
  local reference_name=$4
  local output_suffix=$5

  # Run Minimap2 and save the SAM file in the SAM directory
  minimap2 -a "$reference_genome" "$query_sequence" > "$output_dir/sam/${reference_name}_minimap2${output_suffix}.sam"

  # Filter the best matching result for each gene
  filter_best_match "$output_dir/sam/${reference_name}_minimap2${output_suffix}.sam" "$output_dir/sam/${reference_name}_filtered_minimap2${output_suffix}.sam"

  # Find and realign soft-clipped regions
  extract_and_realign_soft_clipped "$output_dir/sam/${reference_name}_filtered_minimap2${output_suffix}.sam" "$reference_genome" "$output_dir" "$reference_name" "$output_suffix"

  # Convert SAM file into a BED file and save it in the BED directory
  bedtools bamtobed -i "$output_dir/sam/${reference_name}_filtered_minimap2${output_suffix}.sam" > "$output_dir/bed/${reference_name}_minimap2${output_suffix}.bed"

  # Create genome length information file and save it in the BED directory
  cat "$reference_genome" | awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' | sed 's/>//; s/ / /' > "$output_dir/bed/${reference_name}_genome_size${output_suffix}.txt"
}


# Iterate through files in the reference directory
for reference_genome in "$reference_dir"/*.fasta; do
  reference_name=$(basename "$reference_genome" | cut -d '.' -f 1)
  if [ -e "$output_dir/sam/${reference_name}_minimap2.sam" ]; then
    echo "File $reference_name already exists. Skipping."
  else
    echo "Processing $reference_name"

    # First run of the process
    process_alignment "$query_sequence" "$reference_genome" "$output_dir" "$reference_name" ""

    # Second run of the process using soft-clipped sequences from the first run
    process_alignment "$output_dir/temp/${reference_name}_soft_clipped.fasta" "$reference_genome" "$output_dir" "$reference_name" "_second"

    # Third run of the process using soft-clipped sequences from the second run
    process_alignment "$output_dir/temp/${reference_name}_soft_clipped_second.fasta" "$reference_genome" "$output_dir" "$reference_name" "_third"

    # Combine the final SAM files
    # Write header from the first SAM file to the combined SAM file
    grep "^@" "$output_dir/sam/${reference_name}_filtered_minimap2_third.sam" > "$output_dir/sam/${reference_name}_combined.sam"
    # Append the alignment records from the first, second, and third runs to the combined SAM file
    grep -v "^@" "$output_dir/sam/${reference_name}_filtered_minimap2.sam" >> "$output_dir/sam/${reference_name}_combined.sam"
    grep -v "^@" "$output_dir/sam/${reference_name}_filtered_minimap2_second.sam" >> "$output_dir/sam/${reference_name}_combined.sam"
    grep -v "^@" "$output_dir/sam/${reference_name}_filtered_minimap2_third.sam" >> "$output_dir/sam/${reference_name}_combined.sam"
  fi
done
