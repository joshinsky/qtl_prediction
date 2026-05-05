#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

from transformers import AutoTokenizer, AutoModel, AutoConfig
import torch
import sys
import pandas as pd
import numpy as np
import os
import psutil
from huggingface_hub import hf_hub_download
from transformers.models.bert.configuration_bert import BertConfig
import h5py

# get user input
if len(sys.argv) < 3:
    raise SystemExit("Usage: ./DNABERT2_embeddings.py <dataset_name> <window_len>")
dataset = sys.argv[1]
window_len = sys.argv[2]


# get embeddings
def get_embeds(sequences, output_file, batch_size, max_length, device):
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


# load raw string sequences
source_filename = dataset
seq_df = pd.read_csv(source_filename, compression='gzip', usecols=[f'variant_window_{window_len}_ref', f'variant_window_{window_len}_alt'], sep='\t', low_memory=False)
ref_sequences = list(seq_df[f'variant_window_{window_len}_ref'])
alt_sequences = list(seq_df[f'variant_window_{window_len}_alt'])
print(f"\nProcessing {len(seq_df)} sequences...")

# load the Tokenizer and Config
print("\nLoading model...")
tokenizer = AutoTokenizer.from_pretrained("quietflamingo/dnabert2-no-flashattention", trust_remote_code=True)
config = BertConfig.from_pretrained("quietflamingo/dnabert2-no-flashattention")
config.pad_token_id = tokenizer.pad_token_id

# build model with pre-trained weights
model = AutoModel.from_pretrained("quietflamingo/dnabert2-no-flashattention", trust_remote_code=True)
model_path = hf_hub_download(repo_id="quietflamingo/dnabert2-no-flashattention", filename="pytorch_model.bin")
device = torch.device("cpu")
model = model.to(device)
model.eval()

# get maximum token length
approx_seq_len = (int(window_len) * 2) + 50
max_length = int(approx_seq_len * 0.25) + 50    # account for Byte-pair encoding
if max_length > 512:
    max_length = 512
print(f"Using max token length of {max_length} for {window_len}bp window...")

# get embeddings
print("\nCreating embeddings...")
batch_size = 250
ref_output_file = f"results/output/dataset_prep/ref_{window_len}_embeddings_DNABERT2.h5"
alt_output_file = f"results/output/dataset_prep/alt_{window_len}_embeddings_DNABERT2.h5"
get_embeds(ref_sequences, ref_output_file, batch_size, max_length, device)
get_embeds(alt_sequences, alt_output_file, batch_size, max_length, device)
