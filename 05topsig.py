#!/home/ctools/opt/anaconda-2025-12-2/bin/python3

import numpy as np
import pandas as pd
import sys

source_filename = sys.argv[1]
destination = sys.argv[2]
df = pd.read_csv(source_filename, compression='gzip', sep='\t')

# get most significant variants per gene
print(f"Extracting the most significantly associated variants for each gene from {source_filename}...")
df_sorted = df.sort_values('pvalue', ascending=True)
most_sig = df_sorted.drop_duplicates(subset='gene_id', keep='first')
most_sig.to_csv(destination, sep='\t', compression='gzip', index=False)
print(f"Analysis finished with {len(most_sig)} resulting variants selected!\nStored at {destination}")