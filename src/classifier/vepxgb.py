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

# convert h5 embeds into numpy arrays
with h5py.File("results/output/classifier/train_embeddings.h5", 'r') as h5_in:
	X_train = h5_in['embeddings'][...]
with h5py.File("results/output/classifier/validation_embeddings.h5", 'r') as h5_in:
	X_val = h5_in['embeddings'][...]

# load significance data
y_train = pd.read_csv("results/output/classifier/train_dataset.tsv.gz", usecols=['sig_ge', 'sig_iu'] , compression="gzip", sep='\t', low_memory=False)
y_val = pd.read_csv("results/output/classifier/validation_dataset.tsv.gz", usecols=['sig_ge', 'sig_iu'], compression="gzip", sep='\t', low_memory=False)


####################
## PCA on X_train ##
####################

# scale X_train
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_val = sc.transform(X_val)

# Fit PCA on training data
n_components = 300
pca = PCA(n_components=n_components, random_state=42)
X_train_pca = pca.fit_transform(X_train)
X_val_pca = pca.transform(X_val)

explained_variance = pca.explained_variance_ratio_.sum()
print(f"Total explained variance by {n_components} PCs: {explained_variance:.4f}")


#############################
## Plot Explained Variance ##
#############################

# Get the explained variance ratio for each PC
explained_variance_ratio = pca.explained_variance_ratio_

# Calculate the cumulative explained variance
cumulative_variance = np.cumsum(explained_variance_ratio)

# Create the plot
plt.figure(figsize=(10, 6))

# Bar chart for individual explained variance
plt.bar(range(1, n_components + 1), explained_variance_ratio, alpha=0.6, color='b',
        label='Individual Explained Variance')

# Step plot (line) for cumulative explained variance
plt.step(range(1, n_components + 1), cumulative_variance, where='mid', color='r',
         label='Cumulative Explained Variance')

# Aesthetics
plt.title('PCA Explained Variance', fontsize=14)
plt.ylabel('Explained Variance Ratio')
plt.xlabel('Principal Component Index')
plt.xticks(np.arange(0, n_components + 1, step=10)) # Adjust step size if needed
plt.legend(loc='best')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

# Save the plot
output_plot_path = "results/output/classifier/pca_explained_variance.png"
plt.savefig(output_plot_path, dpi=300)
print(f"PCA variance plot saved to {output_plot_path}")