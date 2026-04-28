#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import pandas as pd
import numpy as np
import sys

# load dataset
df = pd.read_csv("results/output/dataset_prep/deduplicated_dataset.tsv.gz", sep='\t', compression='gzip')

def inject_variant(row):
    # variant ID is formatted as 'chr_pos_ref_alt'
    var_parts = row['variant'].split('_')
    ref_allele = var_parts[2].upper()
    alt_allele = var_parts[3].upper()
    
    # Get the wild-type sequence
    ref_seq = row['variant_window'].upper()
    
    # ensure the reference allele actually matches the sequence at index 100
    expected_ref_in_seq = ref_seq[100:100+len(ref_allele)]
    if expected_ref_in_seq != ref_allele:
        print("something doesn't match!")
        sys.exit(1)
        # return pd.NA
        
    # splice in the alternative allele
    # sequence before (first 100bp) + ALT allele + sequence after
    alt_seq = ref_seq[:100] + alt_allele + ref_seq[100+len(ref_allele):]
    
    return alt_seq

print("Injecting ALT alleles to create variant sequences...")
df['variant_window_alt'] = df.apply(inject_variant, axis=1)

# Drop any rows where the sanity check failed (if necessary)
df = df.dropna(subset=['variant_window_alt'])

print("Saving updated dataset...")
# Save the new dataset
df.to_csv("results/output/dataset_prep/final_full_dataset.tsv.gz", sep='\t', index=False, compression='gzip')
print("Done!")
