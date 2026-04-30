#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import sys
import joblib
from sklearn.metrics import roc_auc_score


####################
## get user input ##
####################

try: 
	chosen_classifier = sys.argv[1]
	pca_components = sys.argv[2]
except IndexError:
	print("not enough input arguments. Usage:\npython3 final_evaluation.py <chosen_classifier> <pca_components>")
	sys.exit(1)

print("Loading saved model...")
try:
	model = joblib.load(f"results/output/classifier/trained_{chosen_classifier}_model.joblib")
except FileNotFoundError:
	print("could not find specified model file. Double-check that it exists.")
	sys.exit(1)

if pca_components != 'skipPCA':
	try:	
		sc = joblib.load("results/output/classifier/fitted_scaler.joblib")
		pca = joblib.load(f"results/output/classifier/fitted_pca_{pca_components}_comp.joblib")
	except FileNotFoundError:
		print("could not find all specified files. Double-check that they exist.")
		sys.exit(1)


##################
# load test set ##
##################

def get_available_memory():
	try:
		with open('/proc/meminfo', 'r') as f:
			for line in f:
				if 'MemAvailable' in line:
					return int(line.split()[1]) / (1024**2)
	except Exception:
		return 0.0
	return 0.0

def load_embedding(filename, test_indices):
	with h5py.File(filename, 'r') as embeds_in:
		dataset = embeds_in['embeddings']
		available_RAM = get_available_memory()
		print(f'available system memory: {available_RAM:.1f} GB')

		if available_RAM >= 8.0:
			print("that's enough to load all embeddings quickly!")
			embeddings = dataset[...]
			X_test = embeddings[test_indices]
			del embeddings
		else:
			print("that's not enough and I'll have to use fancy indexing to get embeds.\nThis will be a little slower.")
			X_test = dataset[test_indices]

	return X_test

# load full dataset
test_df = pd.read_csv("results/output/classifier/test_dataset.tsv.gz", compression="gzip", sep='\t', low_memory=False)
test_indices = sorted(test_df['source_row'].tolist())

# load embeddings
ref_embeds_path = "results/output/dataset_prep/ref_embeddings_DNABERT2.h5"
alt_embeds_path = "results/output/dataset_prep/alt_embeddings_DNABERT2.h5"
ref_emb_test = load_embedding(ref_embeds_path, test_indices)
alt_emb_test = load_embedding(alt_embeds_path, test_indices)

print("Calculating Delta Embeddings (ALT - REF)...")
X_test = alt_emb_test - ref_emb_test

# create X and y
test_mapping = np.searchsorted(test_indices, test_df['source_row'].tolist())
X_test = X_test[test_mapping][...]
y_test = test_df[['sig_ge', 'sig_iu']]
y_test_np = y_test.values


##############
# apply PCA ##
##############

if pca_components != 'skipPCA':
	X_test_scaled = sc.transform(X_test)
	X_test_pca = pca.transform(X_test_scaled)
else:
	X_test_pca = X_test


#############
# evaluate ##
#############

print("Running predictions on held-out test set...")
test_probs = model.predict_proba(X_test_pca)
test_prob_sig_ge = test_probs[0][:, 1]
test_prob_sig_iu = test_probs[1][:, 1]

# calculate and print AUC
auc_test_ge = roc_auc_score(y_test_np[:, 0], test_prob_sig_ge)
auc_test_iu = roc_auc_score(y_test_np[:, 1], test_prob_sig_iu)
print(f"Final Test AUC - Gene Expression: {auc_test_ge:.4f}")
print(f"Final Test AUC - Isoform Usage: {auc_test_iu:.4f}")
