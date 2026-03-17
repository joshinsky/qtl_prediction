#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import numpy as np
import pandas as pd
import sys
import time
start_time = time.time()

source_filename = sys.argv[1]
destination = sys.argv[2]
df = pd.read_csv(source_filename, compression='gzip', sep='\t', low_memory=False)

# get most significant variants per gene
print(f"Extracting most significant variants per gene from {source_filename}...")
# Force pvalue column to be numeric, turning any strings into NaNs
df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
df_sorted = df.sort_values('pvalue', ascending=True)
most_sig = df_sorted.drop_duplicates(subset='gene_id', keep='first')
most_sig.to_csv(destination, sep='\t', compression='gzip', index=False)
total_time = time.time() - start_time
print(f'finished after {total_time/60:.2f} minutes!')
print(f"selected {len(most_sig)} variants.\nStored at {destination}\n")