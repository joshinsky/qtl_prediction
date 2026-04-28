#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM python3

import pandas as pd
import h5py
import os


def store_split(split_df, tsv_output_path, h5_output_path, source_h5_path):
    
	print(f"\nProcessing split saving to: {tsv_output_path}")

	###########################
	##### DIAGNOSIS ###########
	###########################
	print(f"split_df rows: {len(split_df)}")
	print(f"split_df index min/max: {split_df.index.min()} / {split_df.index.max()}")
	print(f"split_df head indices: {split_df.index[:5].tolist()}")
	print(f"split_df tail indices: {split_df.index[-5:].tolist()}")


	split_df.to_csv(tsv_output_path, compression="gzip", sep='\t', header=True, index=False)

	# store embeddings using seq-frame indeces
	indices = split_df.index.tolist()

	###########################
	##### DIAGNOSIS ###########
	###########################
	print(f"indices length: {len(indices)}")
	print(f"indices min/max: {min(indices)} / {max(indices)}")
	print(f"indices head: {indices[:5]}")
	print(f"indices tail: {indices[-5:]}")

	print(f"Extracting {len(indices)} corresponding embeddings...")
	with h5py.File(source_h5_path, 'r') as h5_in, h5py.File(h5_output_path, 'w') as h5_out:
		dataset_in = h5_in['embeddings']

		###########################
		##### DIAGNOSIS ###########
		###########################
		print(f"source H5 shape: {dataset_in.shape}")
		if len(indices) > 0:
			print(f"first source embedding row shape: {dataset_in[indices[0]].shape}")

		dataset_out = h5_out.create_dataset("embeddings", shape=(len(indices), 768), dtype='float32', compression="gzip")
		dataset_out[:] = dataset_in[indices]

	###########################
	##### DIAGNOSIS ###########
	###########################
	with h5py.File(h5_output_path, "r") as h5_check:
		print(f"written H5 shape: {h5_check['embeddings'].shape}")

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
output_test_seqs = "results/output/classifier/test_dataset.tsv.gz"
output_train_val_seqs = "results/output/classifier/train_and_validation_dataset.tsv.gz"
input_h5 = "results/output/dataset_prep/embeddings_DNABERT2.h5"
output_test_h5 = "results/output/classifier/test_embeddings.h5"
output_train_val_h5 = "results/output/classifier/train_and_validation_dataset.h5"

###########################
##### DIAGNOSIS ###########
###########################
print(f"test_df rows: {len(test_df)}")
print(f"train_val_df rows: {len(train_val_df)}")

if not os.path.exists(output_test_h5) or not os.path.exists(output_test_seqs):
	store_split(test_df, output_test_seqs, output_test_h5, input_h5)
if not os.path.exists(output_train_val_h5) or not os.path.exists(output_train_val_seqs):
	store_split(train_val_df, output_train_val_seqs, output_train_val_h5, input_h5)

# reset index so val/train splits index correctly into train_and_validation_dataset.h5
train_val_df = train_val_df.reset_index(drop=True)

# look at the data distribution in the set for training + validation
print("\ndistribution data for train+val subset:")
get_chr_distribution(train_val_df)


print("\nI will hold out chr2, chr5, and chr16 for a validation set-size of 20,01%")
val_mask = (train_val_df['chromosome'] == 'chr2') | (train_val_df['chromosome'] == 'chr5') | (train_val_df['chromosome'] == 'chr16')
val_df = train_val_df[val_mask]
train_df = train_val_df[~val_mask]

# # define all directories and perform train-validation split
output_val_seqs = "results/output/classifier/validation_dataset.tsv.gz"
output_train_seqs = "results/output/classifier/train_dataset.tsv.gz"
output_val_h5 = "results/output/classifier/validation_embeddings.h5"
output_train_h5 = "results/output/classifier/train_embeddings.h5"

###########################
##### DIAGNOSIS ###########
###########################
print(f"val_df rows: {len(val_df)}")
print(f"train_df rows: {len(train_df)}")


if not os.path.exists(output_val_h5) or not os.path.exists(output_val_seqs):
	store_split(val_df, output_val_seqs, output_val_h5, output_train_val_h5)
if not os.path.exists(output_train_h5) or not os.path.exists(output_train_seqs):
	store_split(train_df, output_train_seqs, output_train_h5, output_train_val_h5)

print("\nAll splits successfully processed and saved!")

