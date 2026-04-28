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
from transformers.models.bert.configuration_bert import BertConfig
import h5py


start_time = time.time()

if len(sys.argv) < 2:
    raise SystemExit("Usage: ./DNABERT2_embeddings.py <dataset_name>")
dataset = sys.argv[1]

source_filename = f"results/output/dataset_prep/final_full_dataset.tsv.gz"
seq_df = pd.read_csv(source_filename, compression='gzip', usecols=['variant_window_alt'], sep='\t', low_memory=False)
sequences = list(seq_df['variant_window_alt'])
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
max_length = 120    # based on max sequence length of 453 bp in dataset (Byte-Pair Encoding: enc_length = 0.25*seq_length)
batch_size = 2000
all_mean_embeddings = []
output_file = f"results/output/dataset_prep/alt_embeddings_DNABERT2.h5"

# get embeddings
with h5py.File(output_file, 'w') as h5f:
    emb_dataset = h5f.create_dataset(
        "embeddings", 
        shape=(0, 768), 
        maxshape=(None, 768), 
        dtype='float32',
        compression="gzip")

    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch_seqs = sequences[i:i+batch_size]
        
            encoded_inputs = tokenizer(batch_seqs, return_tensors="pt", padding="max_length", max_length=max_length)
            input_ids = encoded_inputs["input_ids"].to(device)
            attention_mask = encoded_inputs["attention_mask"].to(device)
        
            # compute embeddings
            outputs = model(input_ids, attention_mask=attention_mask)
            embeddings = outputs[0]
        
            # Expand mask for math operations
            mask_expanded = attention_mask.unsqueeze(-1)
            mean_batch_embeddings = ( (embeddings * mask_expanded).sum(dim=1) / mask_expanded.sum(dim=1).clamp(min=1) )
        
            # Convert to numpy CPU array
            batch_np = mean_batch_embeddings.cpu().numpy()

            # make space in emb_dataset for new entry & then add it
            current_size = emb_dataset.shape[0]
            emb_dataset.resize(current_size + batch_np.shape[0], axis=0)
            emb_dataset[current_size:] = batch_np

            if current_size%20 == 0:
                print(f"Processed {current_size + batch_np.shape[0]} / {len(sequences)}", end='\r')

print(f"Saved embeddings to {output_file}")