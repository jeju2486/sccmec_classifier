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

def within_distance(position, mec_pos, genome_size, distance=20000):
    if position is None:
        return False
    
    distance_forward = abs(position - mec_pos)
    distance_backward = genome_size - distance_forward
    
    return min(distance_forward, distance_backward) <= distance

def determine_mec_type(mec_pos, mec_strand, mecRI_pos, mecI_pos, IS431_positions, IS1272_pos, blaZ_pos, genome_size, distance=20000):
    if mec_pos is None:
        return "no_sccmec"

    if mec_strand == '+':
        IS431_up_pos = next((pos for pos in IS431_positions if pos < mec_pos), None)
        IS431_down_pos = next((pos for pos in IS431_positions if pos > mec_pos), None)
    else:
        IS431_up_pos = next((pos for pos in IS431_positions if pos > mec_pos), None)
        IS431_down_pos = next((pos for pos in IS431_positions if pos < mec_pos), None)
    
    print(f"mec_pos: {mec_pos}, mec_strand: {mec_strand}, mecRI_pos: {mecRI_pos}, mecI_pos: {mecI_pos}, IS431_up_pos: {IS431_up_pos}, IS431_down_pos: {IS431_down_pos}, IS1272_pos: {IS1272_pos}, blaZ_pos: {blaZ_pos}")

    if mec_strand == '-':
      if blaZ_pos is not None and blaZ_pos < mec_pos and within_distance(blaZ_pos, mec_pos, genome_size, distance):
        print(f"Type E detected for mec_pos {mec_pos}: blaZ_pos {blaZ_pos}")
        return "typeE"
    
      elif (within_distance(IS431_up_pos, mec_pos, genome_size, distance) and within_distance(IS431_down_pos, mec_pos, genome_size, distance)):
        print(f"Type C detected for mec_pos {mec_pos}: IS431_up_pos {IS431_up_pos}, IS431_down_pos {IS431_down_pos}")
        return "typeC"
    
      elif IS1272_pos is not None and IS1272_pos > mec_pos and within_distance(IS1272_pos, mec_pos, genome_size, distance):
        print(f"Type B detected for mec_pos {mec_pos}: IS1272_pos {IS1272_pos}")
        return "typeB"

      elif (within_distance(mecRI_pos, mec_pos, genome_size, distance) or
            within_distance(mecI_pos, mec_pos, genome_size, distance) and
            within_distance(IS431_down_pos, mec_pos, genome_size, distance)):
            print(f"Type A detected for mec_pos {mec_pos}: mecRI_pos {mecRI_pos}, mecI_pos {mecI_pos}, IS431_down_pos {IS431_down_pos}")
            return "typeA"
      else:
        return "unknown"

    elif mec_strand == '+':
      if blaZ_pos is not None and blaZ_pos > mec_pos and within_distance(blaZ_pos, mec_pos, genome_size, distance):
        print(f"Type E detected for mec_pos {mec_pos}: blaZ_pos {blaZ_pos}")
        return "typeE"
    
      elif (within_distance(IS431_up_pos, mec_pos, genome_size, distance) and within_distance(IS431_down_pos, mec_pos, genome_size, distance)):
        print(f"Type C detected for mec_pos {mec_pos}: IS431_up_pos {IS431_up_pos}, IS431_down_pos {IS431_down_pos}")
        return "typeC"
    
      elif IS1272_pos is not None and IS1272_pos < mec_pos and within_distance(IS1272_pos, mec_pos, genome_size, distance):
        print(f"Type B detected for mec_pos {mec_pos}: IS1272_pos {IS1272_pos}")
        return "typeB"
      
      elif (within_distance(mecRI_pos, mec_pos, genome_size, distance) or
            within_distance(mecI_pos, mec_pos, genome_size, distance) and
            within_distance(IS431_down_pos, mec_pos, genome_size, distance)):
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
        blaZ_pos = None

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
                elif gene_name == "blaZ":
                    blaZ_pos = position
        
        mec_pos = mecA_pos if mecA_pos is not None else mecC_pos
        mec_strand = mecA_strand if mecA_strand is not None else mecC_strand

        mec_type = determine_mec_type(mec_pos, mec_strand, mecRI_pos, mecI_pos, IS431_positions, IS1272_pos, blaZ_pos, genome_size)

        results.append([file_prefix + "_scaffold", mec_type, mec_pos, mec_strand])  # Ensuring consistent prefix

    return results

def main():
    bed_dir = './bed'
    sam_dir = './sam'
    output_file = 'mec_info.tsv'

    bed_results = process_bed_files(bed_dir)
    df_bed = pd.DataFrame(bed_results, columns=["file_prefix", "mec_type", "mec_pos", "mec_strand"])
    df_bed.to_csv(output_file, sep='\t', index=False, header=["name", "mec_type", "mec_pos", "mec_strand"])

    print("All files are saved")

if __name__ == "__main__":
    main()
