#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM python3

import pandas as pd
import h5py


def store_split(split_df, tsv_output_path, h5_output_path, source_h5_path):
    
	print(f"\nProcessing split saving to: {tsv_output_path}")

	# store sequence-data df
	split_df.to_csv(tsv_output_path, compression="gzip", sep='\t', header=True, index=False)

	# store embeddings using seq-frame indeces
	indices = sorted(split_df.index.tolist())
	print(f"Extracting {len(indices)} corresponding embeddings...")
	with h5py.File(source_h5_path, 'r') as h5_in, h5py.File(h5_output_path, 'w') as h5_out:

		dataset_in = h5_in['embeddings']
		dataset_out = h5_out.create_dataset(
			"embeddings", 
			shape=(len(indices), 768), 
			dtype='float32',
			compression="gzip")
		dataset_out[:] = dataset_in[indices]

	print(f"Successfully saved embeddings to {h5_output_path}")


def get_chr_distribution(df):
	# get chromosome distribution across dataset
	abs_counts = df['chromosome'].value_counts()
	rel_freqs = df['chromosome'].value_counts(normalize=True)
	summary_df = pd.concat([abs_counts, rel_freqs], axis=1, keys=['Absolute Count', 'Relative Freq'])
	print(f"total variants: {sum(abs_counts)}")
	print("variant distribution across chromosomes:")
	print(summary_df)


#################
## run program ##
#################

# load sequence dataset
df_full = pd.read_csv("results/output/dataset_prep/full_dataset_deduplicated.tsv", sep='\t', low_memory=False)

# look at the data distribution across chromosomes
print("\ndistribution data for full dataset")
get_chr_distribution(df_full)

# chr1, chr19 and chr8 roughly add up to 0.2, so let's hold those out
print("\nI will hold out chr1, chr19, and chr8 for a testset-size of 19.97%")
test_mask = (df_full['chromosome'] == 'chr1') | (df_full['chromosome'] == 'chr19') | (df_full['chromosome'] == 'chr8')
test_df = df_full[test_mask]
train_val_df = df_full[~test_mask]

# define all directories and perform hold-out
# output_test_seqs = "results/output/classifier/test_dataset.tsv.gz"
# output_train_val_seqs = "results/output/classifier/train_and_validation_dataset.tsv.gz"
# input_h5 = "results/output/dataset_prep/embeddings_DNABERT2.h5"
# output_test_h5 = "results/output/classifier/test_embeddings.h5"
# output_train_val_h5 = "results/output/classifier/train_and_validation_dataset.h5"

# store_split(test_df, output_test_seqs, output_test_h5, input_h5)
# store_split(train_val_df, output_train_val_seqs, output_train_val_h5, input_h5)


# look at the data distribution in the set for training + validation
print("\ndistribution data for train+val subset:")
get_chr_distribution(train_val_df)


# print("\nI will hold out chr1, chr19, and chr8 for a testset-size of 19.97%")
# val_mask = (df['chromosome'] == 'chr1') | (df['chromosome'] == 'chr19') | (df['chromosome'] == 'chr8')
# val_df = df[val_mask]
# train_df = df[~val_mask]

# # define all directories and perform train-validation split
# output_val_seqs = "results/output/classifier/validation_dataset.tsv.gz"
# output_train_seqs = "results/output/classifier/train_dataset.tsv.gz"
# output_val_h5 = "results/output/classifier/validation_embeddings.h5"
# output_train_h5 = "results/output/classifier/train_embeddings.h5"

# store_split(val_df, output_val_seqs, output_val_h5, output_train_val_h5)
# store_split(train_df, output_train_seqs, output_train_h5, output_train_val_h5)

print("\nAll splits successfully processed and saved!")

