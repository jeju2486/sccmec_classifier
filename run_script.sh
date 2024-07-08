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

conda activate ./sccmec_classifier_env

bash run_distinguish_type.sh -q "$query_sequence" -r "$reference_dir" -o "$output_dir" -s "$save_temp"

echo "process stat file..."

output_file1="$output_dir/search_stat.tsv"

if [ -e "$output_file1" ]; then
    echo "Warning: File $output_file1 already exists. overwriting."
    rm "$output_file1"
fi

# Call the Python script to process the SAM files and combine results into a single TSV file
python process_alignments.py "$output_dir/sam" "$query_sequence" "$output_file1"

if [ -e "$output_file1" ]; then
    echo "File saved in $output_file1"
else
    echo "Running Failed. Check the log"
fi

echo "process finding the best matches..."

output_file2="$output_dir/best_result.tsv"

if [ -e "$output_file2" ]; then
    echo "Warning: File $output_file2 already exists. overwriting."
    rm "$output_file2"
fi

python best_select.py "$output_file1" "$output_file2"

if [ -e "$output_file2" ]; then
    echo "File saved in $output_file2"
else
    echo "Running Failed. Check the log"
fi
