#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import pandas as pd

# load ge and iu datasets
df_ge = pd.read_csv("results/output/dataset_prep/ge_dataset.tsv.gz", compression="gzip", sep='\t', low_memory=False)
df_iu = pd.read_csv("results/output/dataset_prep/iu_dataset.tsv.gz", compression="gzip", sep='\t', low_memory=False)

# find out if a variant was significant at least once in the dataset
sig_ge = df_ge.groupby('variant')['significant'].max().rename('sig_ge')
sig_iu = df_iu.groupby('variant')['significant'].max().rename('sig_iu')
df_sig_labels = pd.concat([sig_ge, sig_iu], axis=1).fillna(0).astype(int).reset_index()

# combine dfs and drop "non_significant" column as it is redundant
df_combined = pd.concat([df_ge, df_iu], ignore_index=True)
if "non_significant" in df_combined.columns:
    df_combined = df_combined.drop(columns=["non_significant"])

# remove duplicate variants and 
df_combined = df_combined.sort_values(by='significant', ascending=False)
df_meta = df_combined.drop_duplicates(subset=['variant'], keep='first')

# drop 'significant' column and replace it with new sig_labels
df_meta = df_meta.drop(columns=['significant'])
df_final = pd.merge(df_meta, df_sig_labels, on='variant', how='inner')

# store deduplicated df
df_final.to_csv("results/output/dataset_prep/full_dataset_deduplicated.tsv.gz", compression="gzip", sep='\t', header=True, index=False)