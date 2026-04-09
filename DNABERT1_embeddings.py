#!/usr/bin/env python3

import os
import sys
import time
import torch
import pandas as pd
import numpy as np
from transformers import BertTokenizer, BertModel

start_time = time.time()

if len(sys.argv) < 3:
    raise SystemExit("Usage: ./DNABERT_embeddings.py <dataset_name> <dnabert_model_path> [k]")

dataset = sys.argv[1]
model_path = sys.argv[2]
k = int(sys.argv[3]) if len(sys.argv) > 3 else 6

source_filename = f"results/{dataset}/{dataset}_dataset.tsv.gz"
seq_df = pd.read_csv(
    source_filename,
    compression="gzip",
    usecols=["variant_window"],
    sep="\t",
    low_memory=False
)
sequences = seq_df["variant_window"].dropna().astype(str).tolist()

print(f"\nProcessing {len(sequences)} sequences from {source_filename}...")
print(f"Using DNABERT model at: {model_path}")
print(f"Using k-mer size: {k}")

def seq_to_kmers(seq, k=6):
    seq = seq.upper()
    if len(seq) < k:
        return seq
    return " ".join(seq[i:i+k] for i in range(len(seq) - k + 1))

print("\nConverting sequences to k-mers...")
kmer_sequences = [seq_to_kmers(seq, k=k) for seq in sequences]

print("\nLoading DNABERT tokenizer/model...")
tokenizer = BertTokenizer.from_pretrained(model_path, do_lower_case=False)
model = BertModel.from_pretrained(model_path)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)
model.eval()

print(f"Running on device: {device}")

batch_size = 32 if torch.cuda.is_available() else 8
max_length = 512
all_mean_embeddings = []

print("\nCreating embeddings...")

with torch.no_grad():
    for i in range(0, len(kmer_sequences), batch_size):
        batch_seqs = kmer_sequences[i:i + batch_size]

        encoded = tokenizer(
            batch_seqs,
            return_tensors="pt",
            padding="max_length",
            truncation=True,
            max_length=max_length
        )

        input_ids = encoded["input_ids"].to(device)
        attention_mask = encoded["attention_mask"].to(device)

        outputs = model(input_ids=input_ids, attention_mask=attention_mask)
        last_hidden = outputs.last_hidden_state

        mask_expanded = attention_mask.unsqueeze(-1)
        mean_batch_embeddings = (
            (last_hidden * mask_expanded).sum(dim=1) /
            mask_expanded.sum(dim=1).clamp(min=1)
        )

        all_mean_embeddings.append(mean_batch_embeddings.cpu().numpy())

        batch_num = i // batch_size + 1
        total_batches = (len(kmer_sequences) + batch_size - 1) // batch_size
        print(f"Processed batch {batch_num}/{total_batches}")

mean_sequence_embeddings = np.vstack(all_mean_embeddings)

output_file = f"results/{dataset}/{dataset}_embeddings_DNABERT.npz"
np.savez_compressed(output_file, embeddings=mean_sequence_embeddings)

print(f"\nFinal embeddings shape: {mean_sequence_embeddings.shape}")
print(f"Saved embeddings to {output_file}")
print(f"Total runtime: {(time.time() - start_time)/60:.2f} min")
