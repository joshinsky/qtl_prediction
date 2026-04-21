#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import pandas as pd
import h5py

# load sequence dataset
df = pd.read_csv("results/output/dataset_prep/full_dataset_deduplicated.tsv", sep='\t', low_memory=False)

# get chromosome distribution across dataset
abs_counts = df['chromosome'].value_counts()
rel_freqs = df['chromosome'].value_counts(normalize=True)
summary_df = pd.concat([abs_counts, rel_freqs], axis=1, keys=['Absolute Count', 'Relative Freq'])
print(f"total variants: {sum(abs_counts)}")
print("variant distribution across chromosomes:")
print(summary_df)

# chr1, chr19 and chr8 roughly add up to 0.2, so let's filter those out
print("\nI will hold out chr1, chr19, and chr8 for a testset-size of 19.97%")
mask = (df['chromosome'] == 'chr1') | (df['chromosome'] == 'chr19') | (df['chromosome'] == 'chr8')
test_df = df[mask]

# store hold-out testset
output_seqs = "results/output/dataset_prep/holdout_dataset.tsv.gz"
test_df.to_csv(output_seqs, compression="gzip", sep='\t', header=True, index=False)

# get corresponding embeddings
print("\nExtracting corresponding embeddings...")

# use indeces from testset
test_indices = sorted(test_df.index.tolist())
input_h5 = "results/output/dataset_prep/embeddings_DNABERT2.h5"
output_h5 = "results/output/dataset_prep/holdout_embeddings.h5"

with h5py.File(input_h5, 'r') as h5_in, h5py.File(output_h5, 'w') as h5_out:

	# create dataset of testsize x embedding length
	dataset_in = h5_in['embeddings']
	dataset_out = h5_out.create_dataset(
        "embeddings", 
        shape=(len(test_indices), 768), 
        dtype='float32',
        compression="gzip"
    )

	# extract correct embeddings
    dataset_out[:] = dataset_in[test_indices]

print(f"Successfully saved {len(test_indices)} embeddings to {output_h5}")




