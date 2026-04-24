#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

###########################
# load train and val set ##
###########################

with h5py.File("results/output/classifier/train_embeddings.h5", 'r') as h5_in:
	X_train = h5_in['embeddings']
with h5py.File("results/output/classifier/validation_embeddings.h5", 'r') as h5_in:
	X_val = h5_in['embeddings']

y_train = pd.read_csv("results/output/classifier/train_dataset.tsv.gz", usecols=['sig_ge', 'sig_iu'] , compression="gzip", sep='\t', low_memory=False)
y_val = pd.read_csv("results/output/classifier/validation_dataset.tsv.gz", usecols=['sig_ge', 'sig_iu'], compression="gzip", sep='\t', low_memory=False)


####################
## PCA on X_train ##
####################

# Convert HDF5 datasets to NumPy arrays
X_train = X_train[...]
X_val = X_val[...]

# scale X_train
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_val = sc.fit_transform(X_val)

# Fit PCA on training data
n_components = 100
pca = PCA(n_components=n_components, random_state=42)
X_train_pca = pca.fit_transform(X_train)
X_val_pca = pca.transform(X_val)

explained_variance = pca.explained_variance_ratio_.sum()
print(f"Total explained variance by {n_components} PCs: {explained_variance:.4f}")


