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

print("\nMake seqs")
source_filename = "results/seqs/Alasoo_1_ge_dataset.tsv.gz"
seq_df = pd.read_csv(source_filename, compression='gzip', usecols=['variant_window'], sep='\t', low_memory=False)
sequences = list(seq_df['variant_window'])
seq_end_time = time.time() 
print(f'{len(sequences)} seqs made in {(seq_end_time - start_time) / 60:.2f} min!')

# Import the tokenizer and the model
print("\nStart loading")

# 1. Load the Tokenizer and Config
tokenizer = AutoTokenizer.from_pretrained("quietflamingo/dnabert2-fixed", trust_remote_code=True)
config = AutoConfig.from_pretrained("quietflamingo/dnabert2-fixed", trust_remote_code=True)
config.pad_token_id = tokenizer.pad_token_id

# 2. Build the model architecture from scratch (Bypasses the "meta" device bug!)
model = AutoModel.from_config(config, trust_remote_code=True)

# 3. Download and inject the pre-trained weights 
model_path = hf_hub_download(repo_id="quietflamingo/dnabert2-fixed", filename="pytorch_model.bin")

# Load the .bin file using standard torch
state_dict = torch.load(model_path, map_location="cpu", weights_only=True)
model.load_state_dict(state_dict, strict=False)

load_end_time = time.time()
print(f'loaded in {(load_end_time - seq_end_time) / 60:.2f} min!')

print("\nGet embeds")
max_length = 64
batch_size = 5
all_mean_embeddings = []
print(f"Processing {len(sequences)} sequences in batches of {batch_size}...")

# Initialize the process monitor
process = psutil.Process(os.getpid())
mem_before = process.memory_info().rss / (1024 ** 3)
print(f"Memory before embeddings: {mem_before:.2f} GB")

# get embeddings
with torch.no_grad():
    for i in range(0, len(sequences), batch_size):
        # Process sequences in small batches
        batch_seqs = sequences[i:i+batch_size]
        
        encoded_inputs = tokenizer(batch_seqs, return_tensors="pt", padding="max_length", max_length=max_length)
        tokens_ids = encoded_inputs["input_ids"]
        attention_mask = encoded_inputs["attention_mask"]
        
        # Compute embeddings for this specific batch
        # FIX: Removed the invalid encoder_attention_mask kwarg
        torch_outs = model(
            tokens_ids,
            attention_mask=attention_mask,
            output_hidden_states=True
        )
        
        embeddings = torch_outs[0].detach()
        
        # Expand mask for math operations
        mask_expanded = torch.unsqueeze(attention_mask, dim=-1)
        
        # Compute mean embeddings per sequence in this batch
        mean_batch_embeddings = torch.sum(mask_expanded * embeddings, axis=-2) / torch.sum(mask_expanded, axis=1)
        
        # Store the batch result as a numpy array
        all_mean_embeddings.append(mean_batch_embeddings.cpu().numpy())
        
        # Track memory usage per batch
        current_mem = process.memory_info().rss / (1024 ** 3)
        print(f"Processed batch {(i//batch_size) + 1} / {(len(sequences)//batch_size) + 1} | RAM: {current_mem:.2f} GB")

# Combine all batch arrays into one final array
mean_sequence_embeddings = np.vstack(all_mean_embeddings)

mem_after = process.memory_info().rss / (1024 ** 3)
print(f"\nFinal memory usage: {mem_after:.2f} GB (Net increase: {mem_after - mem_before:.2f} GB)")

print(f"Final embeddings shape: {mean_sequence_embeddings.shape}")
emb_end_time = time.time()
print(f'embeds made in {(emb_end_time - load_end_time) / 60:.2f} min!')
print(f'Total script time: {(emb_end_time - start_time) / 60:.2f} min!')

# Create an output path
output_file = "results/seqs/Alasoo_1_embeddings_DNABERT2.npz"

# Save as a compressed numpy array
np.savez_compressed(output_file, embeddings=mean_sequence_embeddings)
print(f"Saved embeddings to {output_file}")