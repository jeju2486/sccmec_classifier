import pandas as pd

def load_type_info(type_info_file):
    type_info_df = pd.read_csv(type_info_file, sep='\t')
    type_info_df['key'] = type_info_df['mec_type'] + '_' + type_info_df['ccr_type']
    type_info_dict = type_info_df.set_index('key')['sccmec_type'].to_dict()
    return type_info_dict

def add_type_column(filtered_genome_file, type_info_file, output_file):
    # Load the type info data
    type_info_dict = load_type_info(type_info_file)

    # Load the filtered genome data
    filtered_genome_df = pd.read_csv(filtered_genome_file, sep='\t')

    # Create the key column in the filtered genome data
    filtered_genome_df['key'] = filtered_genome_df['mec_type'] + '_' + filtered_genome_df['ccr_type']

    # Map the sccmec_type from the type_info data to the filtered genome data
    filtered_genome_df['type'] = filtered_genome_df['key'].map(type_info_dict).fillna('no_sccmec')

    # Drop the key column as it is no longer needed
    filtered_genome_df.drop('key', axis=1, inplace=True)

    # Save the updated dataframe to a new TSV file
    filtered_genome_df.to_csv(output_file, sep='\t', index=False)

def main():
    filtered_genome_file = 'filtered_genome_type_info.tsv'
    type_info_file = 'type_info.tsv'
    output_file = 'genome_type_info.tsv'

    add_type_column(filtered_genome_file, type_info_file, output_file)
    print("The file has been updated and saved.")

if __name__ == "__main__":
    main()
