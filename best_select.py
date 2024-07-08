import sys
import pandas as pd

def read_data(input_file):
    return pd.read_csv(input_file, sep='\t')

def get_best_mec(row):
    mecA_ratio = float(row['mecA_match_ratio']) if pd.notna(row['mecA_match_ratio']) else 0.0
    mecC_ratio = float(row['mecC_match_ratio']) if pd.notna(row['mecC_match_ratio']) else 0.0

    if mecA_ratio > 0.90 or mecC_ratio > 0.90:
        if mecA_ratio >= mecC_ratio:
            return 'mecA', mecA_ratio
        else:
            return 'mecC', mecC_ratio
    else:
        return 'no_sccmec', 0.0

def get_best_class(row, mec_contigs, mec='no_sccmec'):
    class_columns = [
        'IS1272-dmecR1_classB_typeI_match_ratio',
        'mecI-mecR1_classA_match_ratio',
        'IS1272-dmecR1_classB_typeIV_match_ratio',
        'IS431-dmecR1_classC2_typeV_match_ratio',
        'IS1272-dmecR1_classB_typeVI_match_ratio',
        'IS431-dmecR1_classC1_typeVII_match_ratio',
        'IS431-dmecR1-classC2_typeIX_match_ratio',
        'IS431-dmecR1_classC1_typeX_match_ratio',
        'blaZ-mecA-mecR1-mecI_classE_typeXI_match_ratio',
        'IS431-dmecR1_classC2_typeXII_match_ratio',
        'IS1182-dmecI-mecR1_classA.3_typeII_match_ratio',
        'dmecI-IS1182-dmecI-mecR1_classA.4_typeII_match_ratio'
    ]
    best_class = None
    best_ratio = 0.0
    second_best_class = None
    second_best_ratio = 0.0

    for col in class_columns:
        if pd.notna(row[col]):
            ratio = float(row[col])
            if ratio > best_ratio:
                second_best_class, second_best_ratio = best_class, best_ratio
                best_class, best_ratio = col, ratio
            elif ratio > second_best_ratio:
                second_best_class, second_best_ratio = col, ratio
    
    best_contig_list = []
    second_best_contig_list = []

    if best_class is not None:
        best_contig_list = row[f'{best_class.split("_match_ratio")[0]}_contig'].split(',')
    if second_best_class is not None:
        second_best_contig_list = row[f'{second_best_class.split("_match_ratio")[0]}_contig'].split(',')

    if best_ratio - second_best_ratio < 0.1:
        if any(contig in mec_contigs for contig in best_contig_list):
            print(f"Best contig in mec_contigs: {best_contig_list}")
            return best_class.split("_")[1], best_ratio
        elif any(contig in mec_contigs for contig in second_best_contig_list):
            print(f"Second best contig in mec_contigs: {second_best_contig_list}")
            return second_best_class.split("_")[1], second_best_ratio

    if best_class is not None and mec != 'no_sccmec':
        return best_class.split("_")[1], best_ratio
    else:
        return 'no_mec_class', 0.0

def get_best_ccr(row, mec_contigs, mec = 'no_sccmec'):
    
    
    ccrA_types = ['ccrA1', 'ccrA2', 'ccrA3', 'ccrA4']
    ccrB_types = ['ccrB1', 'ccrB2', 'ccrB3', 'ccrB4', 'ccrB6']
    ccrC_types = ['ccrC1', 'ccrC2']
    
    best_ccrA = None
    best_ccrA_ratio = 0.0
    best_ccrA_contig = None

    best_ccrB = None
    best_ccrB_ratio = 0.0
    best_ccrB_contig = None

    best_ccrC = None
    best_ccrC_ratio = 0.0
    best_ccrC_contig = None

    best_ratio=0.0

    for ccr in ccrA_types:
        ccr_ratio=row[f'{ccr}_match_ratio']
        ccr_contig=row[f'{ccr}_contig']
        
        if ccr_ratio >= best_ratio:
            best_ccrA, best_ccrA_ratio, best_ccrA_contig = ccr, ccr_ratio, ccr_contig
            best_ratio = ccr_ratio
    
    best_ratio=0.0
    
    for ccr in ccrB_types:
        ccr_ratio=row[f'{ccr}_match_ratio']
        ccr_contig=row[f'{ccr}_contig']
        
        if ccr_ratio >= best_ratio:
            best_ccrB, best_ccrB_ratio, best_ccrB_contig = ccr, ccr_ratio, ccr_contig
            best_ratio = ccr_ratio
    
    best_ratio=0.0
    
    for ccr in ccrC_types:
        ccr_ratio=row[f'{ccr}_match_ratio']
        ccr_contig=row[f'{ccr}_contig']
        
        if ccr_ratio >= best_ratio:
            best_ccrC, best_ccrC_ratio, best_ccrC_contig = ccr, ccr_ratio, ccr_contig
            best_ratio = ccr_ratio
        
    combined_ccr=''
    
    if mec != 'no_sccmec':
        if abs(best_ccrC_ratio - best_ccrA_ratio) < 0.1 or abs(best_ccrC_ratio - best_ccrB_ratio) < 0.1:
            if best_ccrC_contig in mec_contigs and best_ccrA_contig not in mec_contigs:
                return best_ccrC, best_ccrC_ratio
            else:
                combined_ccr= str(best_ccrA)+str(best_ccrB)
                return combined_ccr, (best_ccrA_ratio+best_ccrB_ratio)/2
        else:
            if best_ccrA_ratio == max(best_ccrA_ratio,best_ccrB_ratio,best_ccrC_ratio) or best_ccrB_ratio == max(best_ccrA_ratio,best_ccrB_ratio,best_ccrC_ratio):
                combined_ccr= str(best_ccrA)+str(best_ccrB)
                return combined_ccr, (best_ccrA_ratio+best_ccrB_ratio)/2
            else:
                return best_ccrC, best_ccrC_ratio     
    else:
        return 'no_ccr', 0.0

def process_data(df):
    result = []
    for index, row in df.iterrows():
        mec_contigs = {row['mecA_contig'], row['mecC_contig']}
        best_mec, best_mec_ratio = get_best_mec(row)
        best_class, best_class_ratio = get_best_class(row, mec_contigs,best_mec)
        best_ccr, best_ccr_ratio = get_best_ccr(row, mec_contigs,best_mec)
        result.append({
            'reference_name': row['sample_id'],
            'best_mec': best_mec,
            'best_mec_ratio': best_mec_ratio,
            'best_mec_class': best_class,
            'best_class_ratio': best_class_ratio,
            'best_ccr': best_ccr,
            'best_ccr_ratio':  best_ccr_ratio
        })
    return pd.DataFrame(result)

def save_results(df, output_file):
    df.to_csv(output_file, sep='\t', index=False)

def main(input_file, output_file):
    df = read_data(input_file)
    processed_df = process_data(df)
    save_results(processed_df, output_file)

if __name__ == "__main__":
    input_file = sys.argv[1] 
    output_file = sys.argv[2] 
    main(input_file, output_file)
