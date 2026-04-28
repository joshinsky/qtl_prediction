#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from flaml import AutoML
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc
from sklearn.multioutput import MultiOutputClassifier


###########################
# load train and val set ##
###########################

# convert h5 embeds into numpy arrays
with h5py.File("results/output/classifier/train_embeddings.h5", 'r') as h5_in:
	X_train = h5_in['embeddings'][...]
with h5py.File("results/output/classifier/validation_embeddings.h5", 'r') as h5_in:
	X_val = h5_in['embeddings'][...]

# load significance data and convert to numpy
y_train = pd.read_csv("results/output/classifier/train_dataset.tsv.gz", usecols=['sig_ge', 'sig_iu'] , compression="gzip", sep='\t', low_memory=False)
y_val = pd.read_csv("results/output/classifier/validation_dataset.tsv.gz", usecols=['sig_ge', 'sig_iu'], compression="gzip", sep='\t', low_memory=False)
y_train_np = y_train.values
y_val_np = y_val.values

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



######################
## Train classifier ##
######################

# load user-specified classifier choice
chosen_classifier = sys.argv[1]

# setup single-label classifier
automl = AutoML()
automl.fit(
	task='classification',
	estimator_list=[chosen_classifier],
	time_budget=30,
	metric='roc_auc',
	verbose=0
	)

# wrap single-label classifier into multi-output classifer
model = MultiOutputClassifier(automl)
model.fit(X_train_pca, y_train_np)

# predict and extract probabilities for each class
train_probs = model.predict_proba(X_train_pca)
val_probs = model.predict_proba(X_val_pca)
train_prob_sig_ge = train_probs[0][:, 1]
val_prob_sig_ge = val_probs[0][:, 1]
train_prob_sig_iu = train_probs[1][:, 1]
val_prob_sig_iu = val_probs[1][:, 1]

# Calculate ROC AUC
auc_train_ge = roc_auc_score(y_train_np[:, 0], train_prob_sig_ge)
auc_val_ge = roc_auc_score(y_val_np[:, 0], val_prob_sig_ge)
auc_train_iu = roc_auc_score(y_train_np[:, 1], train_prob_sig_iu)
auc_val_iu = roc_auc_score(y_val_np[:, 1], val_prob_sig_iu)
print(f"ROC AUC for Gene Expression (sig_ge) - Train: {auc_train_ge:.4f}, Val: {auc_val_ge:.4f}")
print(f"ROC AUC for Isoform Usage (sig_iu) - Train: {auc_train_iu:.4f}, Val: {auc_val_iu:.4f}")



##################
## Plot metrics ##
##################

def plot_roc(y_true_np, y_probs, target_names, title, out_path):
	plt.figure(figsize=(8,6))
	colors = ['blue', 'green']

	# plot ROC for each targets
	for i, target in enumerate(target_names):
		prob_significant = y_probs[i][:, 1]

		# get FPR, TPR and AUC
		fpr, tpr, _ = roc_curve(y_true_np[:, i], prob_significant)
		roc_auc = auc(fpr, tpr)

		# plot curve
		plt.plot(fpr, tpr, color=colors[i], lw=2,
			label=f'{target} (AUC = {roc_auc:.3f})'
			)


	# plot baseline
	plt.plot([0,1], [0,1], color='gray', lw=2, linestyle='--', label='baseline model')

	# edit plot look
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('FPR', fontsize=12)
	plt.ylabel('TPR', fontsize=12)
	plt.title(f'ROC Curve - {title}', fontsize=14)
	plt.legend(loc="lower right")
	plt.grid(alpha=0.3)

	# save and show
	plt.tight_layout()
	plt.savefig(out_path, dpi=500)
	plt.close()


# plot ROCs
plot_roc(
	y_true_np=y_train_np,
	y_probs=train_probs,
	target_names=['significant expression-level effect', 'significant isoform-level effect'],
	title='ROC for training',
	out_path='results/figures/roc_curve_train.png'
	)
plot_roc(
	y_true_np=y_val_np,
	y_probs=val_probs,
	target_names=['significant expression-level effect', 'significant isoform-level effect'],
	title='ROC for validation',
	out_path='results/figures/roc_curve_validation.png'
	)








