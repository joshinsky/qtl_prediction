#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from flaml import AutoML, tune
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay, precision_recall_curve, average_precision_score
from sklearn.multioutput import MultiOutputClassifier
import joblib
import argparse


####################
## get user input ##
####################

def parse_arguments():
	parser = argparse.ArgumentParser(description="Run XGBoost on genomic embeddings")
	parser.add_argument("--classifier", type=str, default="xgboost", choices=["xgboost"])
	parser.add_argument("--pca", type=str, default="skip", help="Provide an int (e.g. 300), 'auto' (for 95% variance), or 'skip'")
	parser.add_argument("--window_size", type=str, choices=["20", "100", "1000"], required=True)
	parser.add_argument("--embedding_type", type=str, choices=["alt", "delta"], required=True)
	parser.add_argument("--target_label", type=str, choices=["standard", "single"], help="'standard' for multi-output, 'single' for two separate classifiers", default="standard")
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
	# val_df = pd.read_csv("results/output/classifier/validation_dataset.tsv.gz", compression="gzip", sep='\t', low_memory=False)
	val_df = pd.read_csv("results/output/classifier/test_dataset.tsv.gz", compression="gzip", sep='\t', low_memory=False)		# uncomment for final eval
	train_indices = sorted(train_df['source_row'].tolist())
	val_indices = sorted(val_df['source_row'].tolist())

	# load embedding files based on window size
	ref_embeds_path = f"results/output/dataset_prep/ref_{args.window_size}_embeddings_DNABERT2.h5"
	alt_embeds_path = f"results/output/dataset_prep/alt_{args.window_size}_embeddings_DNABERT2.h5"
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
	# sc = StandardScaler()
	# X_train = sc.fit_transform(X_train)
	sc = joblib.load("results/output/classifier/xgboost_wt-weighted_tgt-standard_pca-auto_win-1000_emb-alt/xgboost_wt-weighted_tgt-standard_pca-auto_win-1000_emb-alt_fitted_scaler.joblib")
	X_train = sc.transform(X_train)
	X_val = sc.transform(X_val)

	# Save the scaler
	# if outfile_name != '-':
	# 	joblib.dump(sc, f"results/output/classifier/{outfile_name}_fitted_scaler.joblib")

	# Fit PCA on training data
	# if n_components == 'auto':
	# 	pca = PCA(n_components=0.95, random_state=42)
	# else:
	# 	pca = PCA(n_components=int(n_components), random_state=42)

	# X_train_pca = pca.fit_transform(X_train)
	pca = joblib.load("results/output/classifier/xgboost_wt-weighted_tgt-standard_pca-auto_win-1000_emb-alt/xgboost_wt-weighted_tgt-standard_pca-auto_win-1000_emb-alt_fitted_pca_auto_comp.joblib")
	X_train_pca = pca.transform(X_train)
	X_val_pca = pca.transform(X_val)

	actual_n_components = pca.n_components_	

	# Save the PCA model
	# if outfile_name != '-':
	# 	joblib.dump(pca, f"results/output/classifier/{outfile_name}_fitted_pca_{n_components}_comp.joblib")

	# Get the explained variance ratio for each PC
	explained_variance = pca.explained_variance_ratio_.sum()
	print(f"Total explained variance by {actual_n_components} PCs: {explained_variance:.4f}")
	explained_variance_ratio = pca.explained_variance_ratio_

	# Calculate the cumulative explained variance
	cumulative_variance = np.cumsum(explained_variance_ratio)

	# Create the plot
	plt.figure(figsize=(10, 6))

	# Bar chart for individual explained variance
	plt.bar(range(1, actual_n_components + 1), explained_variance_ratio, alpha=0.6, color='b',
        	label='Individual Explained Variance')

	# Step plot (line) for cumulative explained variance
	plt.step(range(1, actual_n_components + 1), cumulative_variance, where='mid', color='r',
		label='Cumulative Explained Variance')

	# Aesthetics
	plt.title('PCA Explained Variance', fontsize=14)
	plt.ylabel('Explained Variance Ratio')
	plt.xlabel('Principal Component Index')

	# dynamically adjust step-size based on n_components
	step_size = max(1, actual_n_components // 10)	
	plt.xticks(np.arange(0, actual_n_components + 1, step=step_size))
	
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
			"scale_pos_weight": {"domain": tune.uniform(lower=1.0, upper=average_scale_weight * 1.2)},
			# "scale_pos_weight": {"domain": average_scale_weight},
			"max_depth": {"domain": tune.randint(lower=3, upper=10)},
			"min_child_weight": {"domain": tune.randint(lower=1, upper=10)},
			"subsample": {"domain": tune.uniform(lower=0.5, upper=0.9)}, 
			"colsample_bytree": {"domain": tune.uniform(lower=0.5, upper=0.9)},
			"reg_alpha": {"domain": tune.loguniform(lower=1e-4, upper=10.0)},
			"reg_lambda": {"domain": tune.loguniform(lower=1e-4, upper=10.0)}
		}}

	# # custom hyperparameter space for xgboost			# uncomment this and early stopping for fighting overfitting
	# custom_search_space = {
	# 	args.classifier: {
	# 		"scale_pos_weight": {"domain": tune.uniform(lower=1.0, upper=average_scale_weight * 1.2)},
	# 		# "scale_pos_weight": {"domain": average_scale_weight},
	# 		"max_depth": {"domain": tune.randint(lower=3, upper=6)},
	# 		"min_child_weight": {"domain": tune.randint(lower=1, upper=50)},
	# 		"subsample": {"domain": tune.uniform(lower=0.5, upper=0.9)}, 
	# 		"colsample_bytree": {"domain": tune.uniform(lower=0.5, upper=0.9)},
	# 		"reg_alpha": {"domain": tune.loguniform(lower=1e-4, upper=10.0)},
	# 		"reg_lambda": {"domain": tune.loguniform(lower=1e-4, upper=10.0)}
	# 	}}		

	automl = AutoML(
		task='classification',
		estimator_list=[args.classifier],
		time_budget=3600,
		metric='roc_auc',
		custom_hp=custom_search_space,
		verbose=0
	)

	print("fit model...")
	if args.target_label == 'single':
		print("Training single-output model for sig_ge...")
		model_ge = AutoML(task='classification', estimator_list=[args.classifier], time_budget=3000, metric='roc_auc', custom_hp=custom_search_space, verbose=0)
		model_ge.fit(X_train_pca, y_train_np[:, 0])		# Column sig_ge
		# model_ge.fit(X_train_pca, y_train_np[:, 0], eval_set=[(X_val_pca, y_val_np[:, 0])], early_stopping_rounds=50)		# uncomment for early stopping
		print("Best hyperparameters for GE:", model_ge.best_config)

		print("Training single-output model for sig_iu...")
		model_iu = AutoML(task='classification', estimator_list=[args.classifier], time_budget=3600, metric='roc_auc', custom_hp=custom_search_space, verbose=0)
		model_iu.fit(X_train_pca, y_train_np[:, 1]) 	# Column sig_iu
		# model_iu.fit(X_train_pca, y_train_np[:, 1], eval_set=[(X_val_pca, y_val_np[:, 1])], early_stopping_rounds=50)		# uncomment for early stopping
		print("Best hyperparameters for IU:", model_iu.best_config)

		# store both models in a dict
		model = {'ge': model_ge, 'iu': model_iu}

	else:
		# wrap single-label classifier into multi-output classifer
		model = MultiOutputClassifier(automl)
		model.fit(X_train_pca, y_train_np)
		print("Best hyperparameters for target 0 (GE):", model.estimators_[0].best_config)
		print("Best hyperparameters for target 1 (IU):", model.estimators_[1].best_config)
		

	return model



##################
## Plot metrics ##
##################

def plot_ROC(y_true_np, y_probs, target_names, title, out_path):
	plt.figure(figsize=(8,6))
	colors = ['blue', 'green']

	# plot ROC for each targets
	for i, target in enumerate(target_names):

		prob_significant = y_probs[i][:, 1] 
		y_true_col = y_true_np[:, i]
		
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

def plot_multi_label_ROC(y_true_np, y_probs, title, out_path):
	# Extract probabilities
	prob_ge = y_probs[0][:, 1]
	prob_iu = y_probs[1][:, 1]
	
	true_ge = y_true_np[:, 0]
	true_iu = y_true_np[:, 1]

	# Map true labels to 4 distinct classes
	def map_to_class(ge, iu):
		return ge * 1 + iu * 2

	y_true_mapped = map_to_class(true_ge, true_iu)
	class_names = ['Neither', 'Only GE', 'Only IU', 'Both']

	# Calculate combined probabilities assuming independence
	prob_class_0 = (1 - prob_ge) * (1 - prob_iu) # Neither
	prob_class_1 = prob_ge * (1 - prob_iu)       # Only GE
	prob_class_2 = (1 - prob_ge) * prob_iu       # Only IU
	prob_class_3 = prob_ge * prob_iu             # Both
	
	# Stack into shape (n_samples, 4) and binarize true labels
	y_score_mapped = np.column_stack((prob_class_0, prob_class_1, prob_class_2, prob_class_3))
	y_true_binarized = label_binarize(y_true_mapped, classes=[0, 1, 2, 3])
	
	plt.figure(figsize=(8,6))
	colors = ['gray', 'blue', 'green', 'purple']
	
	fpr_grid = np.linspace(0.0, 1.0, 1000)
	mean_tpr = np.zeros_like(fpr_grid)
	
	for i, color in zip(range(4), colors):
		fpr, tpr, _ = roc_curve(y_true_binarized[:, i], y_score_mapped[:, i])
		roc_auc = auc(fpr, tpr)
		plt.plot(fpr, tpr, color=color, lw=2, label=f'{class_names[i]} (AUC = {roc_auc:.3f})')
		mean_tpr += np.interp(fpr_grid, fpr, tpr)

	mean_tpr /= 4
	macro_auc = auc(fpr_grid, mean_tpr)
	plt.plot(fpr_grid, mean_tpr, color='black', linestyle=':', lw=3, label=f'Macro-average (AUC = {macro_auc:.3f})')

	plt.plot([0, 1], [0, 1], 'k--', lw=2)
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate', fontsize=12)
	plt.ylabel('True Positive Rate', fontsize=12)
	plt.title(f'Multi-class OvR ROC Curve - {title}', fontsize=14)
	plt.legend(loc="lower right")
	plt.grid(alpha=0.3)
	
	roc_out_path = out_path.replace(".png", "_Combined_4x4_ROC.png")
	plt.savefig(roc_out_path, dpi=500)
	plt.close()
	print(f"Saved Combined 4x4 ROC to {roc_out_path}")


def plot_CM(y_true_np, y_probs, target_names, title, out_path, cutoff=0.5):
	for i, target in enumerate(target_names):

		# case: multiple vs one target
		prob_significant = y_probs[i][:, 1]
		y_true_col = y_true_np[:, i]

		# calculate optimal cutoff using Youden's J statistic
		fpr, tpr, thresholds = roc_curve(y_true_col, prob_significant)
		best_idx = np.argmax(tpr - fpr)
		best_cutoff = thresholds[best_idx]
		
		# create confusion matrix object
		y_pred_class = (prob_significant >= best_cutoff).astype(int)
		cm = confusion_matrix(y_true_col, y_pred_class)
		
		# plot confusion matrix
		disp = ConfusionMatrixDisplay(confusion_matrix=cm)
		disp.plot(cmap=plt.cm.Blues)
		plt.title(f"CM ({target}) - {title} (Cutoff: {best_cutoff})")
		
		# save figure
		target_out_path = out_path.replace(".png", f"_{target.replace(' ', '_')}.png")
		plt.savefig(target_out_path, bbox_inches='tight')
		plt.close()
		print(f"Saved Confusion Matrix to {target_out_path}")

def plot_multi_label_CM(y_true_np, y_probs, title, out_path, cutoff=0.5):
	# Extract predictions for GE (index 0) and IU (index 1)
	prob_ge = y_probs[0][:, 1]
	prob_iu = y_probs[1][:, 1]
	
	true_ge = y_true_np[:, 0]
	true_iu = y_true_np[:, 1]

	# Calculate optimal cutoffs independently for GE and IU
	fpr_ge, tpr_ge, thresholds_ge = roc_curve(true_ge, prob_ge)
	best_cutoff_ge = thresholds_ge[np.argmax(tpr_ge - fpr_ge)]

	fpr_iu, tpr_iu, thresholds_iu = roc_curve(true_iu, prob_iu)
	best_cutoff_iu = thresholds_iu[np.argmax(tpr_iu - fpr_iu)]

	# Apply specific cutoffs
	pred_ge = (prob_ge >= best_cutoff_ge).astype(int)
	pred_iu = (prob_iu >= best_cutoff_iu).astype(int)

	# case: defined cutoff
	# pred_ge = (prob_ge >= cutoff).astype(int)
	# pred_iu = (prob_iu >= cutoff).astype(int)

	# Map combinations to 4 distinct classes:
	def map_to_class(ge, iu):
		return ge * 1 + iu * 2

	y_true_mapped = map_to_class(true_ge, true_iu)
	y_pred_mapped = map_to_class(pred_ge, pred_iu)

	# Generate 4x4 Confusion Matrix
	cm = confusion_matrix(y_true_mapped, y_pred_mapped, labels=[0, 1, 2, 3])
	class_names = ['Neither', 'Only GE', 'Only IU', 'Both']
	
	disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=class_names)
	
	plt.figure(figsize=(8,6))
	disp.plot(cmap=plt.cm.Purples, values_format='d')
	plt.title(f"Combined 4x4 CM - {title} (Cutoff: GE: {best_cutoff_ge}, IU: {best_cutoff_iu})")
	plt.xticks(rotation=45)
	
	target_out_path = out_path.replace(".png", "_Combined_4x4.png")
	plt.savefig(target_out_path, bbox_inches='tight')
	plt.close()
	print(f"Saved Combined 4x4 Confusion Matrix to {target_out_path}")

def plot_PR(y_true_np, y_probs, target_names, title, out_path):
	plt.figure(figsize=(8,6))
	colors = ['blue', 'green']

	for i, target in enumerate(target_names):
		prob_significant = y_probs[i][:, 1]
		y_true_col = y_true_np[:, i]
		
		# get Precision, Recall and Average Precision (PR-AUC)
		precision, recall, _ = precision_recall_curve(y_true_col, prob_significant)
		pr_auc = average_precision_score(y_true_col, prob_significant)
		
		# Baseline is the ratio of positive cases
		baseline = y_true_col.sum() / len(y_true_col)

		# plot them
		plt.plot(recall, precision, color=colors[i % len(colors)], lw=2,
			label=f'{target} (PR AUC = {pr_auc:.3f})'
			)
		# plot baseline
		plt.axhline(y=baseline, color=colors[i % len(colors)], lw=1, linestyle='--', alpha=0.5, label=f'{target} Baseline ({baseline:.3f})')

	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('Recall', fontsize=12)
	plt.ylabel('Precision', fontsize=12)
	plt.title(f'Precision-Recall Curve - {title}', fontsize=14)
	plt.legend(loc="upper right")
	plt.grid(alpha=0.3)

	plt.tight_layout()
	plt.savefig(out_path, dpi=500)
	plt.close()
	print(f"successfully created and saved PR curve at:\n{out_path}")

def plot_multi_label_PR(y_true_np, y_probs, target_names, title, out_path):

	# Build combined 4-class probabilities (mirrors plot_multi_label_ROC)
	prob_ge = y_probs[target_names.index("sig_ge")][:, 1]
	prob_iu = y_probs[target_names.index("sig_iu")][:, 1]

	class_probs = {
		"Neither":  (1 - prob_ge) * (1 - prob_iu),
		"Only GE":  prob_ge       * (1 - prob_iu),
		"Only IU":  (1 - prob_ge) * prob_iu,
		"Both":     prob_ge       * prob_iu,
	}

	ge_true = y_true_np[:, target_names.index("sig_ge")]
	iu_true = y_true_np[:, target_names.index("sig_iu")]

	class_true = {
		"Neither":  ((ge_true == 0) & (iu_true == 0)).astype(int),
		"Only GE":  ((ge_true == 1) & (iu_true == 0)).astype(int),
		"Only IU":  ((ge_true == 0) & (iu_true == 1)).astype(int),
		"Both":     ((ge_true == 1) & (iu_true == 1)).astype(int),
	}

	colors = {"Neither": "grey", "Only GE": "blue", "Only IU": "green", "Both": "purple"}

	fig, ax = plt.subplots(figsize=(8, 6))
	ap_scores = {}

	for cls_name, y_score in class_probs.items():
		y_true_cls = class_true[cls_name]
		baseline   = y_true_cls.mean()
		precision, recall, _ = precision_recall_curve(y_true_cls, y_score)
		ap = average_precision_score(y_true_cls, y_score)
		ap_scores[cls_name] = ap

		ax.plot(recall, precision, color=colors[cls_name], lw=2,
			label=f"{cls_name} (AP = {ap:.3f})")
		ax.axhline(baseline, color=colors[cls_name], lw=1, linestyle="--", alpha=0.5,
			label=f"{cls_name} Baseline ({baseline:.3f})")

	macro_ap = np.mean(list(ap_scores.values()))
	ax.axhline(macro_ap, color="black", lw=1.5, linestyle=":",
		label=f"Macro-avg AP ({macro_ap:.3f})")

	ax.set_xlabel("Recall", fontsize=12)
	ax.set_ylabel("Precision", fontsize=12)
	ax.set_title(title, fontsize=13)
	ax.set_xlim([0.0, 1.0])
	ax.set_ylim([0.0, 1.05])
	ax.legend(loc="upper right", fontsize=9)
	ax.grid(True, alpha=0.3)
	plt.tight_layout()
	plt.savefig(out_path, dpi=150)
	plt.close()



def evaluate_results(X_train_pca, X_val_pca, y_train_np, y_val_np, train_df, val_df, model, args, target_names):
	print(f'get predicted probabilities...')

	# 1. Get probabilities for both sets
	if isinstance(model, dict):
		# case: two single classifiers
		train_probs_ge = model['ge'].predict_proba(X_train_pca)
		train_probs_iu = model['iu'].predict_proba(X_train_pca)
		train_probs = [train_probs_ge, train_probs_iu]

		val_probs_ge = model['ge'].predict_proba(X_val_pca)
		val_probs_iu = model['iu'].predict_proba(X_val_pca)
		val_probs = [val_probs_ge, val_probs_iu]
	else:
		# case: multi-output classifier
		train_probs = model.predict_proba(X_train_pca)
		val_probs = model.predict_proba(X_val_pca)

	positions = ['all', 'exonic', 'intronic', 'intergenic', 'intragenic']

	for pos in positions:
		# Setup position masks for BOTH train and val dataframes
		if pos == 'all':
			mask_val = np.ones(len(val_df), dtype=bool)
			mask_train = np.ones(len(train_df), dtype=bool)
			title_suffix = "All"
		elif pos == 'intragenic':
			mask_val = val_df['variant_location'].isin(['exonic', 'intronic'])
			mask_train = train_df['variant_location'].isin(['exonic', 'intronic'])
			title_suffix = "Intragenic"
		else:
			mask_val = val_df['variant_location'] == pos
			mask_train = train_df['variant_location'] == pos
			title_suffix = pos.capitalize()
		
		if not mask_val.any() or not mask_train.any():
			print(f"Warning: No data found for variant_location == {pos}")
			return
			
		# Subset Val data
		y_val_subset = y_val_np[mask_val]
		val_probs_subset = [vp[mask_val] for vp in val_probs] if isinstance(val_probs, list) else val_probs[mask_val]

		# Subset Train data
		y_train_subset = y_train_np[mask_train]
		train_probs_subset = [tp[mask_train] for tp in train_probs] if isinstance(train_probs, list) else train_probs[mask_train]
	
		# Determine file prefixes
		base_out_val = f"results/figures/{args.outfile}_{pos}.png"
		base_out_train = f"results/figures/{args.outfile}_train_{pos}.png"

		# Generate Validation Plots
		if "roc" in args.outfile.lower():
			plot_ROC(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out_val.replace(".png", "_roc.png"))
			plot_multi_label_ROC(y_val_subset, val_probs_subset, f'Validation ({title_suffix})', base_out_val)
		elif "cm" in args.outfile.lower():
			plot_CM(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out_val.replace(".png", "_cm.png"))
			plot_multi_label_CM(y_val_subset, val_probs_subset, f'Validation ({title_suffix})', base_out_val)
		else:
			plot_ROC(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out_val.replace(".png", "_roc.png"))
			plot_CM(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out_val.replace(".png", "_cm.png"))
			plot_multi_label_ROC(y_val_subset, val_probs_subset, f'Validation ({title_suffix})', base_out_val)
			plot_multi_label_CM(y_val_subset, val_probs_subset, f'Validation ({title_suffix})', base_out_val)
			plot_ROC(y_train_subset, train_probs_subset, target_names, f'Train ({title_suffix})', base_out_train.replace(".png", "_roc.png"))
			plot_multi_label_ROC(y_train_subset, train_probs_subset, f'Train ({title_suffix})', base_out_train)
			plot_PR(y_train_subset, train_probs_subset, target_names, f'Train ({title_suffix})', base_out_train.replace(".png", "_PR.png"))
			plot_PR(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out_val.replace(".png", "_PR.png"))
			plot_multi_label_PR(y_train_subset, train_probs_subset, target_names, f'Train ({title_suffix})', base_out_train.replace(".png", "_4x4_PR.png"))
			plot_multi_label_PR(y_val_subset, val_probs_subset, target_names, f'Validation ({title_suffix})', base_out_val.replace(".png", "_4x4_PR.png"))

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
			# Pass 'auto' directly, or convert to int
			pca_param = 'auto' if args.pca.lower() == 'auto' else int(args.pca)
			X_train_pca, X_val_pca = run_pca(X_train, X_val, args.outfile, n_components=pca_param)
		except ValueError:
			print("input number of pca components as int, 'auto', or 'skip'.")
			sys.exit(1)

	# train classifier
	# model = run_classifier(X_train_pca, y_train_np, args)
	model = joblib.load("results/output/classifier/xgboost_wt-weighted_tgt-standard_pca-auto_win-1000_emb-alt/xgboost_wt-weighted_tgt-standard_pca-auto_win-1000_emb-alt.joblib")	# uncomment for final evaluation
	
	# save model
	# if args.outfile != '-':
	# 	model_filename = f"results/output/classifier/{args.outfile}.joblib"
	# 	joblib.dump(model, model_filename)
	# 	print(f"\nModel successfully saved to {model_filename}")

	# evaluate and plot
	target_names = ['sig_ge', 'sig_iu']
	evaluate_results(X_train_pca, X_val_pca, y_train_np, y_val_np, train_df, val_df, model, args, target_names)

	print(f"\nFinished!")


if __name__ == "__main__":
	main()



