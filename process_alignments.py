import os
import sys
import pandas as pd
import pysam

def parse_fasta(file_path):
    gene_lengths = {}
    gene_order = []
    with open(file_path, 'r') as f:
        gene_name = ""
        sequence = ""
        for line in f:
            if line.startswith('>'):
                if gene_name:
                    gene_lengths[gene_name] = len(sequence)
                    gene_order.append(gene_name)
                gene_name = line[1:].strip()
                sequence = ""
            else:
                sequence += line.strip()
        if gene_name:
            gene_lengths[gene_name] = len(sequence)
            gene_order.append(gene_name)
    return gene_lengths, gene_order

def get_match_length(cigartuples):
    return sum(length for code, length in cigartuples if code == 0)  # CIGAR code 0 is match

def parse_sam(file_path):
    matches = {}
    with pysam.AlignmentFile(file_path, "r") as samfile:
        for read in samfile:
            if read.is_unmapped or read.mapping_quality <= 50:
                continue
            query_name = read.query_name
            contig_name = samfile.get_reference_name(read.reference_id)
            match_length = get_match_length(read.cigartuples)
            if query_name not in matches:
                matches[query_name] = {}
            if contig_name not in matches[query_name]:
                matches[query_name][contig_name] = 0
            matches[query_name][contig_name] += match_length
    return matches

def collect_gene_stats(sam_dir):
    gene_stats = {}
    for filename in os.listdir(sam_dir):
        if filename.endswith('_combined.sam'):
            sample_id = filename.split('_')[0]
            file_path = os.path.join(sam_dir, filename)
            matches = parse_sam(file_path)
            gene_stats[sample_id] = matches
    return gene_stats

def format_output(gene_stats, query_lengths, gene_order):
    data = []
    for sample_id, genes in gene_stats.items():
        row = {'sample_id': sample_id}
        for gene in gene_order:
            if gene in genes:
                total_match_length = sum(genes[gene].values())
                contig_names = ",".join(genes[gene].keys())
                gene_total = query_lengths.get(gene, 0)
                match_ratio = f"{(total_match_length / gene_total):.3f}" if gene_total > 0 else "0.000"
                row[f'{gene}_match'] = f"{total_match_length}/{gene_total}"
                row[f'{gene}_match_ratio'] = match_ratio
                row[f'{gene}_contig'] = contig_names
            else:
                row[f'{gene}_match'] = "0/0"
                row[f'{gene}_match_ratio'] = "0.000"
                row[f'{gene}_contig'] = ""
        data.append(row)
    
    df = pd.DataFrame(data)
    return df

def save_to_tsv(df, output_file):
    for column in df.columns:
        if '_match_ratio' in column:
            df[column]=df[column].fillna('0.000')
        elif '_match' in column:
            df[column]=df[column].fillna('0/0')
        elif '_contig' in column:
            df[column]=df[column].fillna('')
    df.to_csv(output_file, sep='\t', index=False)

def main(sam_dir, query_sequence, output_file):
    if not os.path.exists(sam_dir):
        print(f"Error: Directory '{sam_dir}' does not exist.")
        sys.exit(1)
    if not os.path.isfile(query_sequence):
        print(f"Error: Query sequence file '{query_sequence}' does not exist.")
        sys.exit(1)
    
    query_lengths, gene_order = parse_fasta(query_sequence)
    gene_stats = collect_gene_stats(sam_dir)
    df = format_output(gene_stats, query_lengths, gene_order)
    save_to_tsv(df, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_alignments.py <sam_dir> <query_sequence> <output_file>")
        sys.exit(1)
    
    sam_dir = sys.argv[1]
    query_sequence = sys.argv[2]
    output_file = sys.argv[3]
   
    if not output_file:
        print("Error: Output file path cannot be empty.")
        sys.exit(1)
    
    main(sam_dir, query_sequence, output_file)
