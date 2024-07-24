import pandas as pd

# Read the TSV files into DataFrames
df2 = pd.read_csv('type_info.tsv', sep='\t')
df1 = pd.read_csv('best_result.tsv', sep='\t')

# Merge the DataFrames based on the mec_type and ccr_type
merged_df = df1.merge(df2, how='left', left_on=['best_mec_class', 'best_ccr'], right_on=['mec_type', 'ccr_type'])

# Drop the redundant columns from df2
merged_df = merged_df.drop(columns=['mec_type', 'ccr_type'])

# Save the result to a new TSV file
merged_df.to_csv('output.tsv', sep='\t', index=False)

