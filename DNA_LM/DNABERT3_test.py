#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

from transformers import AutoTokenizer, AutoModel
import torch
import time
import sys
import pandas as pd
import numpy as np
import os
import psutil

start_time = time.time()

print("\nMake seqs")
source_filename = "results/seqs/Alasoo_1_ge_dataset.tsv.gz"
seq_df = pd.read_csv(source_filename, compression='gzip', usecols=['variant_window'], sep='\t', low_memory=False)

# NEW: Helper function to convert raw DNA sequence to 3-mer strings with spaces
def seq2kmer(seq, k=3):
    return " ".join([seq[x:x+k] for x in range(len(seq)+1-k)])

# Convert all sequences into the 3-mer format required by DNABERT-1
sequences = [seq2kmer(seq, k=3) for seq in seq_df['variant_window']]

seq_end_time = time.time() 
print(f'{len(sequences)} seqs made in {(seq_end_time - start_time) / 60:.2f} min!')

# Import the tokenizer and the model
print("\nStart loading DNABERT (3-mer)...")
repo_name = "zhihan1996/DNA_bert_3"

# Because this uses standard BERT architecture, we can load it normally!
tokenizer = AutoTokenizer.from_pretrained(repo_name)
model = AutoModel.from_pretrained(repo_name)

load_end_time = time.time()
print(f'loaded in {(load_end_time - seq_end_time) / 60:.2f} min!')

print("\nGet embeds")
# Note: In a 3-mer model, 1 token = 1 nucleotide step. 
# A 250bp sequence yields 248 tokens + 2 special tokens = 250 tokens.
max_length = 256
batch_size = 2000
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
        
        encoded_inputs = tokenizer(batch_seqs, return_tensors="pt", padding="max_length", max_length=max_length, truncation=True)
        tokens_ids = encoded_inputs["input_ids"]
        attention_mask = encoded_inputs["attention_mask"]
        
        # Compute embeddings for this specific batch
        torch_outs = model(
            tokens_ids,
            attention_mask=attention_mask,
            output_hidden_states=True
        )
        
        # Extract the last hidden state (Standard BERT outputs it at index 0)
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
output_file = "results/seqs/Alasoo_1_embeddings_DNABERT3.npz"

# Save as a compressed numpy array
np.savez_compressed(output_file, embeddings=mean_sequence_embeddings)
print(f"Saved embeddings to {output_file}")
