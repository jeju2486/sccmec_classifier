#!/bin/bash

# Function to display usage message
usage() {
  echo "Usage: $0 -q <query_sequence> -r <reference_dir> -o <output_dir> [-s TRUE|FALSE]"
  exit 1
}

# Parse command-line arguments
save_temp="TRUE"
while getopts q:r:o:s: flag
do
    case "${flag}" in
        q) query_sequence=${OPTARG};;
        r) reference_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        s) save_temp=${OPTARG};;
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

# Function to extract genes that start at 1 or end at contig length and save to a new SAM file
extract_genes_for_second_run() {
  local sam_file=$1
  local genome_size_file=$2
  local output_sam_file=$3

  awk -v genome_size_file="$genome_size_file" '
    BEGIN {
      while ((getline < genome_size_file) > 0) {
        contig_lengths[$1] = $2;
      }
    }
    /^@/ { print $0 > "'"$output_sam_file"'"; next }
    !/^@/ {
      contig = $3;
      start = $4;
      match($6, /([0-9]+)M/, arr);
      matched_length = arr[1];
      end = start + matched_length - 1;
      if (start <= 100 || end >= contig_lengths[contig]-100) {
        print $0 > "'"$output_sam_file"'";
      }
    }
  ' "$sam_file"
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

  # Extract genes for second run
  extract_genes_for_second_run "$output_dir/sam/${reference_name}_filtered_minimap2${output_suffix}.sam" "$output_dir/bed/${reference_name}_genome_size.txt" "$output_dir/temp/${reference_name}_for_second_run${output_suffix}.sam"
  
  if [ -s "$output_dir/temp/${reference_name}_for_second_run${output_suffix}.sam" ]; then
    # Find and realign soft-clipped regions
    extract_and_realign_soft_clipped "$output_dir/temp/${reference_name}_for_second_run${output_suffix}.sam" "$reference_genome" "$output_dir" "$reference_name" "$output_suffix"

    # Convert SAM file into a BED file and save it in the BED directory
    bedtools bamtobed -i "$output_dir/sam/${reference_name}_filtered_minimap2${output_suffix}.sam" > "$output_dir/bed/${reference_name}_minimap2${output_suffix}.bed"
  else
    echo "No alignments met criteria for second run"
  fi
}

# Iterate through files in the reference directory
for reference_genome in "$reference_dir"/*.{fasta,fas,fna}; do
  if [ ! -e "$reference_genome" ];then
    echo "Warning: file extension is not correct checking other option..."
    continue
  fi      
  reference_name=$(basename "$reference_genome" | cut -d '.' -f 1)
  if [ -e "$output_dir/sam/${reference_name}_minimap2.sam" ]; then
    echo "File $reference_name already exists. Skipping."
  else
    echo "Processing $reference_name"
    
    # Create genome length information file and save it in the BED directory
    cat "$reference_genome" | awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' | sed 's/>//; s/ / /' > "$output_dir/bed/${reference_name}_genome_size.txt"
    
    # First run of the process
    process_alignment "$query_sequence" "$reference_genome" "$output_dir" "$reference_name" ""

    # Second run of the process using soft-clipped sequences from the first run
    if [ -e "$output_dir/temp/${reference_name}_soft_clipped.fasta" ]; then
      process_alignment "$output_dir/temp/${reference_name}_soft_clipped.fasta" "$reference_genome" "$output_dir" "$reference_name" "_second"
    fi

    # Third run of the process using soft-clipped sequences from the second run
    if [ -e "$output_dir/temp/${reference_name}_soft_clipped_second.fasta" ]; then
      process_alignment "$output_dir/temp/${reference_name}_soft_clipped_second.fasta" "$reference_genome" "$output_dir" "$reference_name" "_third"
    fi

    # Combine the final SAM files
    grep "^@" "$output_dir/sam/${reference_name}_filtered_minimap2.sam" > "$output_dir/sam/${reference_name}_combined.sam"
    grep -v "^@" "$output_dir/sam/${reference_name}_filtered_minimap2.sam" >> "$output_dir/sam/${reference_name}_combined.sam"
    if [ -e "$output_dir/sam/${reference_name}_filtered_minimap2_second.sam" ]; then
        grep -v "^@" "$output_dir/sam/${reference_name}_filtered_minimap2_second.sam" >> "$output_dir/sam/${reference_name}_combined.sam"
        if [ -e "$output_dir/sam/${reference_name}_filtered_minimap2_third.sam" ]; then
            grep -v "^@" "$output_dir/sam/${reference_name}_filtered_minimap2_third.sam" >> "$output_dir/sam/${reference_name}_combined.sam"
        fi
    fi
  fi
done
