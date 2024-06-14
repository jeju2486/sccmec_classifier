import os
import pandas as pd

def load_genome_size(genome_size_file):
    with open(genome_size_file) as f:
        line = f.readline().strip()
        fields = line.split('\t')
        if len(fields) >= 2:
            return int(fields[1])
        else:
            raise ValueError(f"Invalid format in genome size file: {genome_size_file}")

def within_distance(position, mec_pos, genome_size, distance=30000):
    if position is None:
        return False
    
    distance_forward = abs(position - mec_pos)
    distance_backward = genome_size - distance_forward
    
    return min(distance_forward, distance_backward) <= distance

def determine_mec_type(mec_pos, mec_strand, mecRI_pos, mecI_pos, IS431_positions, IS1272_pos, genome_size):
    if mec_pos is None:
        return "no_sccmec"

    if mec_strand == '+':
        IS431_up_pos = next((pos for pos in IS431_positions if pos < mec_pos), None)
        IS431_down_pos = next((pos for pos in IS431_positions if pos > mec_pos), None)
    else:
        IS431_up_pos = next((pos for pos in IS431_positions if pos > mec_pos), None)
        IS431_down_pos = next((pos for pos in IS431_positions if pos < mec_pos), None)
    
    print(f"mec_pos: {mec_pos}, mec_strand: {mec_strand}, mecRI_pos: {mecRI_pos}, mecI_pos: {mecI_pos}, IS431_up_pos: {IS431_up_pos}, IS431_down_pos: {IS431_down_pos}, IS1272_pos: {IS1272_pos}")
    
    if within_distance(IS1272_pos, mec_pos, genome_size):
        print(f"Type B detected for mec_pos {mec_pos}: IS1272_pos {IS1272_pos}")
        return "typeB"
        
    if (within_distance(IS431_up_pos, mec_pos, genome_size) and
        within_distance(IS431_down_pos, mec_pos, genome_size)):
        print(f"Type C detected for mec_pos {mec_pos}: IS431_up_pos {IS431_up_pos}, IS431_down_pos {IS431_down_pos}")
        return "typeC"

    if mec_strand == '-':
        if (within_distance(mecRI_pos, mec_pos, genome_size) or
            within_distance(mecI_pos, mec_pos, genome_size) and
            within_distance(IS431_down_pos, mec_pos, genome_size)):
            print(f"Type A detected for mec_pos {mec_pos}: mecRI_pos {mecRI_pos}, mecI_pos {mecI_pos}, IS431_down_pos {IS431_down_pos}")
            return "typeA"
    else:
        if (within_distance(mecRI_pos, mec_pos, genome_size) or
            within_distance(mecI_pos, mec_pos, genome_size) and
            within_distance(IS431_down_pos, mec_pos, genome_size)):
            print(f"Type A detected for mec_pos {mec_pos}: mecRI_pos {mecRI_pos}, mecI_pos {mecI_pos}, IS431_down_pos {IS431_down_pos}")
            return "typeA"

    return "no_sccmec"

def process_bed_files(bed_dir):
    bed_files = [f for f in os.listdir(bed_dir) if f.endswith('.bed')]
    
    results = []

    for bed_file in bed_files:
        file_prefix = os.path.basename(bed_file).replace('_scaffold.bed', '')
        bed_path = os.path.join(bed_dir, bed_file)
        genome_size_file = os.path.join(bed_dir, f"{file_prefix}_scaffold_genome_size.txt")

        try:
            genome_size = load_genome_size(genome_size_file)
        except ValueError as e:
            print(e)
            continue

        mecA_pos, mecA_strand = None, None
        mecC_pos, mecC_strand = None, None
        mecRI_pos, mecI_pos = None, None
        IS431_positions = []
        IS1272_pos = None

        with open(bed_path) as f:
            for line in f:
                fields = line.strip().split('\t')
                gene_name = fields[3]
                position = int(fields[1])
                strand = fields[5]

                if gene_name == "mecA":
                    mecA_pos = position
                    mecA_strand = strand
                elif gene_name == "mecC":
                    mecC_pos = position
                    mecC_strand = strand
                elif gene_name == "mecRI":
                    mecRI_pos = position
                elif gene_name == "mecI":
                    mecI_pos = position
                elif gene_name == "IS431":
                    IS431_positions.append(position)
                elif gene_name == "IS1272":
                    IS1272_pos = position
        
        mec_pos = mecA_pos if mecA_pos is not None else mecC_pos
        mec_strand = mecA_strand if mecA_strand is not None else mecC_strand

        mec_type = determine_mec_type(mec_pos, mec_strand, mecRI_pos, mecI_pos, IS431_positions, IS1272_pos, genome_size)

        results.append([file_prefix + "_scaffold", mec_type, mec_pos, mec_strand])  # Ensuring consistent prefix

    return results

def process_sam_files(sam_dir):
    sam_files = [f for f in os.listdir(sam_dir) if f.endswith('.sam')]
    
    results = []

    for sam_file in sam_files:
        accession = os.path.basename(sam_file).replace('.sam', '')
        sam_path = os.path.join(sam_dir, sam_file)

        shortest_cigar = None
        shortest_cigar_length = float('inf')
        mecA_position = None

        with open(sam_path) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                fields = line.strip().split('\t')
                qname = fields[0]
                rname = fields[2]
                pos = int(fields[3])
                cigar = fields[5]

                if rname == "*" or cigar == "*":
                    continue
                
                if 'mecA' in qname or 'mecC' in qname:
                    mecA_position = pos
                
                if any(ccr in qname for ccr in ["ccrA1B1", "ccrA2B2", "ccrC1", "ccrC2", "ccrA4", "ccrA1B3", "ccrA1B6", "ccrA3B3", "ccrA4B4"]):
                    ccr_position = pos
                    cigar_length = len(cigar) - sum(c.isdigit() for c in cigar)
                    if mecA_position and within_distance(ccr_position, mecA_position, 5000):
                        if cigar_length < shortest_cigar_length:
                            shortest_cigar_length = cigar_length
                            shortest_cigar = qname

        if shortest_cigar:
            results.append([accession, shortest_cigar])
        else:
            results.append([accession, "no_ccr"])

    # Debug print statements for SAM processing results
    for result in results:
        print(f"SAM result: {result}")

    return results

def combine_results(bed_results, sam_results, output_file):
    df_bed = pd.DataFrame(bed_results, columns=["file_prefix", "mec_type", "mec_pos", "mec_strand"])
    df_sam = pd.DataFrame(sam_results, columns=["accession", "ccr_type"])

    combined = pd.merge(df_bed, df_sam, left_on="file_prefix", right_on="accession", how="inner")
    combined.drop("accession", axis=1, inplace=True)

    # Debug print statement for combined results
    print("Combined DataFrame:")
    print(combined)

    combined.to_csv(output_file, sep='\t', index=False, header=["name", "mec_type", "mec_pos", "mec_strand", "ccr_type"])

def main():
    bed_dir = './bed'
    sam_dir = './sam'
    output_file = 'filtered_genome_type_info.tsv'

    bed_results = process_bed_files(bed_dir)
    sam_results = process_sam_files(sam_dir)
    combine_results(bed_results, sam_results, output_file)
    print("All files are saved")

if __name__ == "__main__":
    main()
