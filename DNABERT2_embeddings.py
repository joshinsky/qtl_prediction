#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

from transformers import AutoTokenizer, AutoModel, AutoConfig
import torch
import time
import sys
import pandas as pd
import numpy as np
import os
import psutil
from huggingface_hub import hf_hub_download

start_time = time.time()

if len(sys.argv) < 2:
    raise SystemExit("Usage: ./DNABERT2_embeddings.py <dataset_name>")
dataset = sys.argv[1]

source_filename = f"results/{dataset}/{dataset}_dataset.tsv.gz"
seq_df = pd.read_csv(source_filename, compression='gzip', usecols=['variant_window'], sep='\t', low_memory=False)
sequences = list(seq_df['variant_window'])
print(f"\nProcessing {len(sequences)} sequences...")

# load the Tokenizer and Config
print("\nLoading model...")
tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
config = BertConfig.from_pretrained("zhihan1996/DNABERT-2-117M")
config.pad_token_id = tokenizer.pad_token_id

# build model with pre-trained weights 
model = AutoModel.from_config(config, trust_remote_code=True)
model_path = hf_hub_download(repo_id="zhihan1996/DNABERT-2-117M", filename="pytorch_model.bin")
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)
model.eval()

print("\nCreating embeddings...")
max_length = 64
batch_size = 1000
all_mean_embeddings = []

# get embeddings
with torch.no_grad():
    for i in range(0, len(sequences), batch_size):
        batch_seqs = sequences[i:i+batch_size]
        
        encoded_inputs = tokenizer(batch_seqs, return_tensors="pt", padding="max_length", max_length=max_length)
        input_ids = encoded_inputs["input_ids"].to(device)
        attention_mask = encoded_inputs["attention_mask"].to(device)
        
        # compute embeddings
        outputs = model(input_ids, attention_mask=attention_mask)
        
        embeddings = outputs.last_hidden_state
        
        # Expand mask for math operations
        mask_expanded = attention_mask.unsqueeze(-1)
        mean_batch_embeddings = ( (embeddings * mask_expanded).sum(dim=1) / mask_expanded.sum(dim=1).clamp(min=1) )
        
        # Store the batch result as a numpy array
        all_mean_embeddings.append(mean_batch_embeddings.cpu().numpy())

# Combine all batch arrays into one final array
mean_sequence_embeddings = np.vstack(all_mean_embeddings)

print(f"Final embeddings shape: {mean_sequence_embeddings.shape}")

# Save as a compressed numpy array
output_file = f"results/{dataset}/{dataset}_embeddings_DNABERT2.npz"
np.savez_compressed(output_file, embeddings=mean_sequence_embeddings)
print(f"Saved embeddings to {output_file}")