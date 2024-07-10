#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Demi Gommers

# Required modules:

"""
"prep_variant_input_fasta_maker.py" <input_file> <output_name>
output_name example: "SMA_test_123"
input_file example: "all_SNVs_pivoted_paraphase_suppl_made_in_R.xlsx"

"""

import pandas as pd
import sys
import os
import re



def save_haplotype_alt_file(df, output):
    """
    Function for splitting and saving dataframe into separate files for further processing
    """
    for sample_hap, group in df.groupby('sample_hap'):
        group.to_csv(f"{output}{sample_hap}.tsv", sep='\t', index=False)
    return df



def pos_ref_alt_file(depth_file, pos_ref_alt_file):
    pos_ref_alt_df = pd.read_csv(pos_ref_alt_file)
    pos_ref_alt_df[['POS', 'REF', 'ALT']] = pos_ref_alt_df['pos_ref_alt'].str.split('_', expand=True)
   
    return pos_ref_alt_df        
            
    
def extract_sample_and_haplotype(input_file_path):
    # Extract the file name from the path
    file_name = os.path.basename(input_file_path)

    # Regular expressions to extract sample name and haplotype
    sample_pattern = r'HG\d+|(SMA\d+)((?:_blood)?(?:_fib_P\d+)*(?:_P\d+)?)'
    haplotype_pattern = r'hap[0-9]+'

    sample_match = re.search(sample_pattern, file_name)
    haplotype_match = re.search(haplotype_pattern, file_name)

    sample_name = sample_match.group(0) if sample_match else None
    haplotype = haplotype_match.group(0) if haplotype_match else None
    
    return sample_name, haplotype    
    

def main():
    input_file = sys.argv[1]
    output_name = sys.argv[2]
    
    all_SNVs_SMA = pd.read_csv(input_file, sep='\t', header=0)
    
    # Perform data manipulation and file saving
    SMA_variants = (all_SNVs_SMA.melt(id_vars=['pos_ref_alt'], var_name='sample_hap')
                  .query("value == 1")
                  .loc[:, ['pos_ref_alt', 'sample_hap']]
                  .groupby('sample_hap')
                  .apply(save_haplotype_alt_file, output=output_name))
                  
  
    print("End of script")

if __name__ == "__main__":
    main()
