import pysam
import os

input_directory = "./sam"
output_file = "ccr_info.tsv"
ccr_types = ["ccrA1B1", "ccrA2B2", "ccrC1", "ccrC2", "ccrA4", "ccrA3B3", "ccrA4B4"]

def calculate_total_length(cigar_string):
    total_length = 0
    number = ''
    for char in cigar_string:
        if char.isdigit():
            number += char
        else:
            if char in 'MIS=X':  # Include relevant operations
                total_length += int(number)
            number = ''
    return total_length

def calculate_aligned_length(cigar_string):
    aligned_length = 0
    number = ''
    for char in cigar_string:
        if char.isdigit():
            number += char
        else:
            if char == 'S':  # Soft-clipped bases
                aligned_length += int(number) * 0.5
            elif char in 'M':  # Include matches, insertions, deletions, etc.
                aligned_length += int(number)
            elif char in 'ID':
                aligned_length += int(number) * (-2)
            number = ''
    return aligned_length

def calculate_highest_ccr_match(input_sam, ccr_types):
    ccr_ratios = {ccr: 0 for ccr in ccr_types}
    mecA_present = 0
    with pysam.AlignmentFile(input_sam, "r") as infile:
        for read in infile:
            if not read.is_unmapped and read.cigarstring:
                aligned_length = calculate_aligned_length(read.cigarstring)
                read_total_length = calculate_total_length(read.cigarstring)
                ratio = aligned_length / read_total_length
                reference_name = read.reference_name
                if read.query_name == 'mecA':
                    mecA_present = 1
                print(f"Read: {read.query_name}, Aligned Length: {aligned_length}, Total Length: {read_total_length}, Ratio: {ratio}, Reference: {reference_name}")
                if read.query_name in ccr_types:
                    ccr_ratios[read.query_name] = ratio

    print(f"CCR Ratios: {ccr_ratios}")
    if mecA_present == 1:
        highest_ccr = max(ccr_ratios, key=ccr_ratios.get)
        highest_ratio = ccr_ratios[highest_ccr]
    else:
        highest_ccr = 'no_ccr'
        highest_ratio = 0
    return highest_ccr, highest_ratio

def process_directory(input_dir, output_file, ccr_types):
    results = []
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".sam"):
            input_sam = os.path.join(input_dir, file_name)
            print(f"Processing {input_sam}...")
            highest_ccr, highest_ratio = calculate_highest_ccr_match(input_sam, ccr_types)
            base_name = os.path.basename(file_name).split('.')[0]
            if highest_ratio < 0.1:
                results.append((base_name, "no_ccr", 0))
            else:
                results.append((base_name, highest_ccr, highest_ratio * 100))
    
    with open(output_file, 'w') as f:
        f.write("name\tCCR_type\tpercentage_match\n")
        for base_name, ccr, percentage in results:
            f.write(f"{base_name}\t{ccr}\t{percentage:.2f}%\n")
    print(f"Results saved to {output_file}")

process_directory(input_directory, output_file, ccr_types)
