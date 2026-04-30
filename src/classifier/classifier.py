#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from flaml import AutoML, tune
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay
from sklearn.multioutput import MultiOutputClassifier
import joblib
import argparse


####################
## get user input ##
####################

def parse_arguments():
	parser = argparse.ArgumentParser(description="Run XGBoost/LightGBM on genomic embeddings")
	parser.add_argument("--classifier", type=str, default="xgboost", choices=["xgboost", "lightgbm"])
	parser.add_argument("--pca", type=str, default="skip", help="Provide an integer like 300, or 'skip'")
	parser.add_argument("--window_size", type=str, choices=["20", "100", "1000"], required=True)
	parser.add_argument("--embedding_type", type=str, choices=["alt", "delta"], required=True)
	parser.add_argument("--gene_position", type=str, choices=["all", "exonic", "intronic", "intergenic"], required=True)
	parser.add_argument("--target_label", type=str, choices=["standard", "both"], help="'standard' for GE+IU, 'both' for overlapping sig", default="standard")
	parser.add_argument("--eval_set", type=str, choices=["val", "test"], default="val")
	parser.add_argument("--class_weighting", type=str, choices=["none", "weighted"], default="none")
	parser.add_argument("--outfile", type=str, required=True)
	return parser.parse_args()


###########################
# load train and val set ##
###########################

def get_available_memory():
	try:
		with open('/proc/meminfo', 'r') as f:
			for line in f:
				if 'MemAvailable' in line:
					return int(line.split()[1]) / (1024**2)
	except Exception:
		return 0.0
	return 0.0

def load_embedding(filename, train_indices, val_indices):
	with h5py.File(filename, 'r') as embeds_in:
		dataset = embeds_in['embeddings']
		available_RAM = get_available_memory()
		print(f'available system memory: {available_RAM:.1f} GB')

		if available_RAM >= 8.0:
			print("that's enough to load all embeddings quickly!")
			embeddings = dataset[...]
			X_train = embeddings[train_indices]
			X_val = embeddings[val_indices]
			del embeddings
		else:
			print("that's not enough and I'll have to use fancy indexing to get embeds.\nThis will be a little slower.")
			X_train = dataset[train_indices]
			X_val = dataset[val_indices]

	return X_train, X_val

def load_data(args):
	# load dataset and extract row indeces
	train_df = pd.read_csv("results/output/classifier/train_dataset.tsv.gz", compression="gzip", sep='\t', low_memory=False)
	val_df = pd.read_csv("results/output/classifier/validation_dataset.tsv.gz", compression="gzip", sep='\t', low_memory=False)
	train_indices = sorted(train_df['source_row'].tolist())
	val_indices = sorted(val_df['source_row'].tolist())

	# load embedding files based on window size
	ref_embeds_path = f"results/output/dataset_prep/ref_embeddings_DNABERT2_{args.window_size}.h5"
	alt_embeds_path = f"results/output/dataset_prep/alt_embeddings_DNABERT2_{args.window_size}.h5"
	alt_emb_train, alt_emb_val = load_embedding(alt_embeds_path, train_indices, val_indices)

	# if delta embedding is desired
	if args.embedding_type == 'delta':
		print("Calculating Delta Embeddings (ALT - REF)...")
		ref_emb_train, ref_emb_val = load_embedding(ref_embeds_path, train_indices, val_indices)
		X_train = alt_emb_train - ref_emb_train
		X_val = alt_emb_val - ref_emb_val
	else:
		X_train = alt_emb_train
		X_val = alt_emb_val

	# create X
	train_mapping = np.searchsorted(train_indices, train_df['source_row'].tolist())
	val_mapping = np.searchsorted(val_indices, val_df['source_row'].tolist())
	X_train = X_train[train_mapping][...]
	X_val = X_val[val_mapping][...]

	# create y
	if args.target_label == 'both':
		train_df['sig_both'] = (train_df['sig_ge'] & train_df['sig_iu']).astype(int)
		val_df['sig_both'] = (val_df['sig_ge'] & val_df['sig_iu']).astype(int)
		y_train = train_df[['sig_both']]
		y_val = val_df[['sig_both']]
	else:
		y_train = train_df[['sig_ge', 'sig_iu']]
		y_val = val_df[['sig_ge', 'sig_iu']]

	# convert to np
	y_train_np = y_train.values
	y_val_np = y_val.values

	return X_train, X_val, y_train_np, y_val_np, train_df, val_df


####################
## PCA on X_train ##
####################

def run_pca(X_train, X_val, outfile_name, n_components=1):

	# scale X_train
	sc = StandardScaler()
	X_train = sc.fit_transform(X_train)
	X_val = sc.transform(X_val)

	# Save the scaler
	if outfile_name != '-':
		joblib.dump(sc, f"results/output/classifier/{outfile_name}_fitted_scaler.joblib")

	# Fit PCA on training data
	pca = PCA(n_components=n_components, random_state=42)
	X_train_pca = pca.fit_transform(X_train)
	X_val_pca = pca.transform(X_val)

	# Save the PCA model
	if outfile_name != '-':
		joblib.dump(pca, f"results/output/classifier/{outfile_name}_fitted_pca_{n_components}_comp.joblib")

    # Get the explained variance ratio for each PC
	explained_variance = pca.explained_variance_ratio_.sum()
	print(f"Total explained variance by {n_components} PCs: {explained_variance:.4f}")
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

	return X_train_pca, X_val_pca


######################
## Train classifier ##
######################

def run_classifier(X_train_pca, y_train_np, args):
	print(f'\nstart training {args.classifier}...')

	# account for class imbalance if desired
	if args.class_weighting == 'none':
		print("Using default weighting (no scaling).")
		average_scale_weight = 1.0
	else:
		# calculate class imbalance
		scale_pos_weights = []
		for i in range(y_train_np.shape[1]):
			pos = y_train_np[:, i].sum()
			neg = len(y_train_np) - pos
			scale_pos_weights.append(neg / pos if pos > 0 else 1)
		
		average_scale_weight = sum(scale_pos_weights) / len(scale_pos_weights)
		print(f"Calculated scale_pos_weight: {average_scale_weight:.2f}")


	# custom hyperparameter space for xgboost
	custom_search_space = {
		args.classifier: {
			"scale_pos_weight": {"domain": average_scale_weight},
			"max_depth": {"domain": tune.randint(lower=2, upper=6)},
			"colsample_bytree": {"domain": tune.uniform(lower=0.4, upper=0.8)},
			"reg_alpha": {"domain": tune.loguniform(lower=0.1, upper=10.0)},
			"reg_lambda": {"domain": tune.loguniform(lower=0.1, upper=10.0)}
		}}

	automl = AutoML(
		task='classification',
		estimator_list=[args.classifier],
		time_budget=300,
		metric='roc_auc',
		custom_hp=custom_search_space,
		verbose=0
	)

	print("fit model...")
	if args.target_label == 'both':
		# train only on one feature
		model = automl
		model.fit(X_train_pca, y_train_np.ravel())
	else:
		# wrap single-label classifier into multi-output classifer
		model = MultiOutputClassifier(automl)
		model.fit(X_train_pca, y_train_np)
		
	return model



##################
## Plot metrics ##
##################

def plot_roc(y_true_np, y_probs, target_names, title, out_path):
	plt.figure(figsize=(8,6))
	colors = ['blue', 'green']

	# plot ROC for each targets
	for i, target in enumerate(target_names):

		# case: multiple vs only one target
		if isinstance(y_probs, list):
			prob_significant = y_probs[i][:, 1] 
		else: 
			rob_significant = y_probs[:, 1]

		# case: multiple vs only one target
		if y_true_np.ndim > 1: 
			y_true_col = y_true_np[:, i] 
		else: 
			y_true_col = y_true_np
		
		# get FPR, TPR and AUC
		fpr, tpr, _ = roc_curve(y_true_col, prob_significant)
		roc_auc = auc(fpr, tpr)

		# plot them
		plt.plot(fpr, tpr, color=colors[i % len(colors)], lw=2,
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

		# save plot
		plt.tight_layout()
		plt.savefig(out_path, dpi=500)
		plt.close()
		print(f"successfully created and saved '{title}' at:\n{out_path}")

def plot_CM(y_true_np, y_probs, target_names, title, out_path, cutoff=0.5):
	for i, target in enumerate(target_names):

		# case: multiple vs one target
		if isinstance(y_probs, list):
			prob_significant = y_probs[i][:, 1]
		else:
			prob_significant = y_probs[:, 1]

		# case: multiple vs one target
		if y_true_np.ndim > 1:
			y_true_col = y_true_np[:, i]
		else:
			y_true_col = y_true_np
		
		# create confusion matrix object
		y_pred_class = (prob_significant >= cutoff).astype(int)
		cm = confusion_matrix(y_true_col, y_pred_class)
		
		# plot confusion matrix
		disp = ConfusionMatrixDisplay(confusion_matrix=cm)
		disp.plot(cmap=plt.cm.Blues)
		plt.title(f"CM ({target}) - {title} (Cutoff: {cutoff})")
		
		# save figure
		target_out_path = out_path.replace(".png", f"_{target.replace(' ', '_')}.png")
		plt.savefig(target_out_path, bbox_inches='tight')
		plt.close()
		print(f"Saved Confusion Matrix to {target_out_path}")

def evaluate_results(X_val_pca, y_val_np, val_df, model, args, target_names):
	print(f'get predicted probabilities...')
	val_probs = model.predict_proba(X_val_pca)
	
	if args.gene_position == 'all':
		mask = np.ones(len(val_df), dtype=bool)
		title_suffix = "All"
	else:
		mask = val_df['position'] == args.gene_position  # Ensure 'Variant_Type' or similar column exists
		title_suffix = args.gene_position.capitalize()
		
		if not mask.any():
			print(f"Warning: No data found for gene_position == {args.gene_position}")
			return
			
	y_val_subset = y_val_np[mask]
	val_probs_subset = [vp[mask] for vp in val_probs] if isinstance(val_probs, list) else val_probs[mask]
	
	# Determine file prefixes
	base_out = f"results/figures/{args.outfile}"
	
	if "roc" in args.outfile.lower():
		plot_ROC(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out)
	elif "cm" in args.outfile.lower():
		plot_CM(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out)
	else:
		plot_ROC(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out.replace(".png", "_roc.png"))
		plot_CM(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out.replace(".png", "_cm.png"))



##################
## Main Program ##
##################

def main():
	args = parse_arguments()
	
	# load Data
	X_train, X_val, y_train_np, y_val_np, train_df, val_df = load_data(args)
	
	# PCA
	if args.pca == 'skip':
		print("skipping pca...")
		X_train_pca, X_val_pca = X_train, X_val
	else:
		try: 
			pca_components = int(args.pca)
			X_train_pca, X_val_pca = run_pca(X_train, X_val, args.outfile, n_components=pca_components)
		except ValueError:
			print("input number of pca components as int or type 'skip'.")
			sys.exit(1)

	# train classifier
	model = run_classifier(X_train_pca, y_train_np, args)
	
	# save model
	if args.outfile != '-':
		model_filename = f"results/output/classifier/{args.outfile}.joblib"
		joblib.dump(model, model_filename)
		print(f"\nModel successfully saved to {model_filename}")

	# evaluate and plot
	target_names = ['Significant Both'] if args.target_label == 'both' else ['sig_ge', 'sig_iu']
	evaluate_results(X_val_pca, y_val_np, val_df, model, args, target_names)

	print(f"\nFinished!")


if __name__ == "__main__":
	main()



