#!/bin/bash

#SBATCH --job-name=test_run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=devel
#SBATCH --time=0-00:03:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/sccmec_classifier/test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/sccmec_classifier/test_output.txt

#python code load
module purge
module load Anaconda3/2023.09-0
source activate $DATA/sccmec_classifier_env

# Set the path to your query sequence
query_sequence="/data/biol-micro-genomics/kell7366/sccmec_classifier/gene_db.fasta"

# Set the directory containing reference sequences and flanking length
reference_dir="/data/biol-micro-genomics/kell7366/sccmec_classifier/Assemblies/"

# Create a directory to store the files
output_dir="/data/biol-micro-genomics/kell7366/sccmec_classifier"

bash run_distinguish_type.sh -q "$query_sequence" -r "$reference_dir" -o "$output_dir"

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
