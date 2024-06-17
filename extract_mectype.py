import pandas as pd
import os

# Paths
type_info_file = 'type_info.tsv'
bed_dir = './bed_direct'
output_file = 'mec_info.tsv'

# Load the type_info.tsv file into a DataFrame
type_info = pd.read_csv(type_info_file, sep='\t')

# Create a dictionary to map sccmec_type to mec_type
type_to_mec = type_info.set_index('sccmec_type')['mec_type'].to_dict()

# Print the type_to_mec dictionary for debugging
print("Type to Mec dictionary:")
print(type_to_mec)

# Create a dictionary to store the mec_type(s) for each isolate
isolate_mec_types = {}
all_isolates = set()

# Process each BED file in the bed_direct directory
for bed_file in os.listdir(bed_dir):
    if bed_file.endswith(".bed"):
        bed_path = os.path.join(bed_dir, bed_file)
        base_name = os.path.basename(bed_file).split('.')[0]  # Extract the base name of the BED file
        all_isolates.add(base_name)  # Track all isolates
        print(f"Processing file: {bed_path}")
        result_bed = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'sccmec_type', 'score', 'strand'])
        
        for index, row in result_bed.iterrows():
            sccmec_type = row['sccmec_type'].strip()  # Strip any extra whitespace
            
            mec_type = type_to_mec.get(sccmec_type)
            print(f"File: {base_name}, SCCmec Type: {sccmec_type}, Mec Type: {mec_type}")
            
            if mec_type is not None:
                if base_name not in isolate_mec_types:
                    isolate_mec_types[base_name] = set()
                isolate_mec_types[base_name].add(mec_type)

# Prepare the final output
output = []
for isolate in all_isolates:
    if isolate in isolate_mec_types:
        mec_types = isolate_mec_types[isolate]
        if len(mec_types) == 1:
            final_mec_type = list(mec_types)[0]
        else:
            final_mec_type = '/'.join(sorted(mec_types))
    else:
        final_mec_type = 'no_mec'
    
    output.append((isolate, final_mec_type))

# Convert output to a DataFrame and save it to a TSV file
output_df = pd.DataFrame(output, columns=['name', 'mec_type'])
output_df.to_csv(output_file, sep='\t', index=False)

print(f"Output saved to {output_file}")
