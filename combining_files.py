import pandas as pd

# Paths to input files
ccr_match_result_file = 'ccr_info.tsv'
isolate_mec_info_file = 'mec_info.tsv'
type_info_file = 'type_info.tsv'
output_file = 'combined_results.tsv'

# Load the TSV files into DataFrames
ccr_match_result = pd.read_csv(ccr_match_result_file, sep='\t')
isolate_mec_info = pd.read_csv(isolate_mec_info_file, sep='\t')
type_info = pd.read_csv(type_info_file, sep='\t')

# Merge the CCR match result with isolate mec info on the "name" column
combined_df = pd.merge(ccr_match_result, isolate_mec_info, on='name', how='outer')

# Merge the combined DataFrame with type_info to get the sccmec_type
# Since the merge is based on both mec_type and ccr_type, we need to ensure these columns exist in combined_df
combined_df = pd.merge(combined_df, type_info[['mec_type', 'ccr_type', 'sccmec_type']], 
                       left_on=['mec_type', 'CCR_type'], right_on=['mec_type', 'ccr_type'], 
                       how='left')

# Drop the extra 'ccr_type' column resulting from the merge
combined_df = combined_df.drop(columns=['ccr_type'])

# Sort the combined DataFrame by the "name" column
combined_df = combined_df.sort_values(by='name')

# Save the combined and sorted DataFrame to a TSV file
combined_df.to_csv(output_file, sep='\t', index=False)

print(f"Combined results saved to {output_file}")
