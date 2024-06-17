import pysam
import os
from pathlib import Path

def calculate_total_length(cigar_string):
    total_length = 0
    number = ''
    for char in cigar_string:
        if char.isdigit():
            number += char
        else:
            if char in 'MIS=X':  # Include soft-clipped bases and other relevant operations
                total_length += int(number)
            number = ''
    return total_length

def filter_and_sort_alignments(input_sam, min_percentage=0.60, min_total_length=5500):
    alignments = []
    with pysam.AlignmentFile(input_sam, "r") as infile:
        for read in infile:
            if not read.is_unmapped:
                aligned_length = read.query_alignment_length
                total_length = calculate_total_length(read.cigarstring)
                if total_length >= min_total_length and total_length != 0:
                    ratio = aligned_length / total_length
                    print(f"Aligned Length: {aligned_length}, Total Length: {total_length}, Ratio: {ratio}")
                    if ratio >= min_percentage:
                        alignments.append((ratio, read))
    alignments.sort(reverse=True, key=lambda x: x[0])  # Sort by ratio in descending order
    return alignments[:2]  # Keep only the top two alignments

def process_directory(input_dir, output_dir, min_percentage=0.60, min_total_length=5500):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".sam"):
            input_sam = os.path.join(input_dir, file_name)
            output_sam = os.path.join(output_dir, file_name)
            print(f"Processing {input_sam}...")
            top_alignments = filter_and_sort_alignments(input_sam, min_percentage, min_total_length)
            with pysam.AlignmentFile(input_sam, "r") as infile, pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:
                for _, alignment in top_alignments:
                    outfile.write(alignment)
            print(f"Filtered SAM written to {output_sam}")

input_directory = "./sam_direct/"
output_directory = "./filtered_sam_direct/"
process_directory(input_directory, output_directory)
