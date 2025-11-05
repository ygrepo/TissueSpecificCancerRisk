import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

def load_and_prepare_data(filepath):
    # Load the data
    data = pd.read_csv(filepath, sep='\t')
    data.rename(columns={data.columns[1]: 'Cancer_Type'}, inplace=True)

    # Debug: Print the first few rows and columns
    print("\nData analyzed: ", filepath)
    print(data.head())
    print(data.columns)

    # Identify all chromosome arm columns (assuming they follow a pattern like 'Xq' or 'Xp')
    chrom_arms = [col for col in data.columns[2:] if 'q' in col or 'p' in col]
    print("\nChr arms analyzed: ", chrom_arms)
    for col in chrom_arms:
        data[col] = pd.to_numeric(data[col], errors='coerce')

    return data, chrom_arms

def perform_permutation_test(data, chrom_arm, num_permutations=10000, random_state=None):
    np.random.seed(random_state)
    p_values = {}

    for cancer_type in data['Cancer_Type'].unique():
        data_cancer = data[data['Cancer_Type'] == cancer_type]
        
        if chrom_arm in data.columns:
            # Handle NA properly by ignoring them in the deletions check
            deletion_cancer = data_cancer[chrom_arm].dropna()
            size_cancer = deletion_cancer.size
            observed_stat = (deletion_cancer <= -1).mean()  # -1 or lower indicates deletion
            
            # Combined pool to sample from
            non_cancer_data = data[data['Cancer_Type'] != cancer_type][chrom_arm].dropna()
            combined_pool = np.concatenate([deletion_cancer, non_cancer_data])
            
            # Generate permutation stats
            permutation_stats = []
            for _ in range(num_permutations):
                np.random.shuffle(combined_pool)
                perm_stat = (combined_pool[:size_cancer] == -1).mean()
                permutation_stats.append(perm_stat)

            # Calculate p-value: Proportion of permutation stats >= observed stat
            p_value = np.mean(np.array(permutation_stats) >= observed_stat)
            p_values[cancer_type] = p_value
        else:
            p_values[cancer_type] = 'NA'  # Chromosome arm not present

    return p_values

def correct_p_values(p_values_dict):
    corrected_p_values_dict = {}
    all_p_values = []

    for chrom_arm in p_values_dict:
        for cancer_type in p_values_dict[chrom_arm]:
            if p_values_dict[chrom_arm][cancer_type] != 'NA':
                all_p_values.append(p_values_dict[chrom_arm][cancer_type])

    # Correct for multiple testing
    reject, corrected_p_values, _, _ = multipletests(all_p_values, method='fdr_bh')
    
    index = 0
    for chrom_arm in p_values_dict:
        corrected_p_values_dict[chrom_arm] = {}
        for cancer_type in p_values_dict[chrom_arm]:
            if p_values_dict[chrom_arm][cancer_type] != 'NA':
                corrected_p_values_dict[chrom_arm][cancer_type] = (p_values_dict[chrom_arm][cancer_type], corrected_p_values[index])
                index += 1
            else:
                corrected_p_values_dict[chrom_arm][cancer_type] = ('NA', 'NA')

    return corrected_p_values_dict

def save_results_to_tsv(results, output_filepath):
    rows = []
    for chrom_arm, cancer_types in results.items():
        for cancer_type, (orig_p_value, corr_p_value) in cancer_types.items():
            rows.append([chrom_arm, cancer_type, orig_p_value, corr_p_value])
    
    results_df = pd.DataFrame(rows, columns=['Chromosome_Arm', 'Cancer_Type', 'Original_P_Value', 'Corrected_P_Value'])
    results_df.to_csv(output_filepath, sep='\t', index=False)

def main(filepath, output_filepath, num_permutations=10000, random_state=None):
    data, chrom_arms = load_and_prepare_data(filepath)
    results = {}
    
    for chrom_arm in chrom_arms:
        p_values = perform_permutation_test(data, chrom_arm, num_permutations, random_state)
        results[chrom_arm] = p_values
    
    corrected_results = correct_p_values(results)
    save_results_to_tsv(corrected_results, output_filepath)


# TCGA full analysis
main('../data/PANCAN_ArmCallsAndAneuploidyScore_092817_filtered.txt',
     'out/TCGA_filtered_results.txt', 
     num_permutations=10000, 
     random_state=42)

# TCGA analysis: no WGD samples
main('../data/PANCAN_ArmCallsAndAneuploidyScore_092817_filtered_noWGD.txt',
     'out/TCGA_filtered_noWGD_results.txt', 
     num_permutations=10000, 
     random_state=42)

# ICGC full analysis
main('../data/broad_values_by_arm.rmcnv.pt_170207_filtered.txt',
     'out/ICGC_filtered_results.txt', 
     num_permutations=10000, 
     random_state=42)

# ICGC analysis: no WGD samples
main('../data/broad_values_by_arm.rmcnv.pt_170207_filtered_noWGD.txt',
     'out/ICGC_filtered_noWGD_results.txt', 
     num_permutations=10000, 
     random_state=42)