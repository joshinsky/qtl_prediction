#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM python3

import pandas as pd
import numpy as np
import h5py
import os

def get_chr_distribution(df, set_name):
	# get chromosome distribution across dataset
	abs_counts = df['chromosome'].value_counts()
	rel_freqs = df['chromosome'].value_counts(normalize=True)
	summary_df = pd.concat([abs_counts, rel_freqs], axis=1, keys=['Absolute Count', 'Relative Freq'])
	print(f"\n--- {set_name} ---")
	print(f"total variants: {sum(abs_counts)}")
	print("distribution across chromosomes:")
	print(summary_df)


###################
## test hold-out ##
###################

# load sequence dataset and store original indeces
df_full = pd.read_csv("results/output/dataset_prep/full_dataset_deduplicated.tsv", sep='\t', low_memory=False)
df_full['source_row'] = df_full.index

# look at the data distribution across chromosomes
get_chr_distribution(df_full, set_name='distribution data for full dataset')

# chr1, chr19 and chr8 roughly add up to 0.2, so let's hold those out
print("\nI will hold out chr1, chr19, and chr8 for a testset-size of 19.97%")
test_mask = (df_full['chromosome'] == 'chr1') | (df_full['chromosome'] == 'chr19') | (df_full['chromosome'] == 'chr8')
test_df = df_full[test_mask]
train_val_df = df_full[~test_mask]


############################
## train-validation split ##
############################

# look at the data distribution in the set for training + validation
get_chr_distribution(train_val_df, set_name='distribution data for train+val subset:')

print("\nI will hold out chr2, chr5, and chr16 for a validation set-size of 20,01%")
val_mask = (train_val_df['chromosome'] == 'chr2') | (train_val_df['chromosome'] == 'chr5') | (train_val_df['chromosome'] == 'chr16')
val_df = train_val_df[val_mask]
train_df = train_val_df[~val_mask]


##################
## store splits ##
##################

os.makedirs("results/output/classifier", exist_ok=True)

test_df.to_csv("results/output/classifier/test_dataset.tsv.gz", compression="gzip", sep='\t', header=True, index=False)
val_df.to_csv("results/output/classifier/validation_dataset.tsv.gz", compression="gzip", sep='\t', header=True, index=False)
train_df.to_csv("results/output/classifier/train_dataset.tsv.gz", compression="gzip", sep='\t', header=True, index=False)

print(f"\nFinal row counts -> Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")
print("\nAll splits successfully processed and saved!")

