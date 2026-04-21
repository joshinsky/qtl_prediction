#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import os
import sys
import gzip
import math
import pandas as pd
import numpy as np
import genomicranges
from genomicranges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
import time
start_time = time.time()

##############################
##   01 - Get user input    ##
## (filenames)              ##
##############################

# function explaining the program usage
def usage():
	print("program usage:")
	print(f"python3 {sys.argv[0]} <filename> <destination>")
	sys.exit(1)

# get user input for filenames
if sys.argv[1] == 'help':
	usage()

# handle filenames
try:
	input_filename = sys.argv[1]
	destination_filename = sys.argv[2]
except IndexError:
	print(f"I am missing a filename.")
	usage()

##############################
## 02 - Load gtf data and   ##
##  get gr objects          ##
##############################

def make_gr_from_df(df):
	meta_df = df.drop(columns=["seqnames", "starts", "ends", "strand"], errors="ignore")
	
	# get positional indeces as numpy arrays
	starts_arr = pd.to_numeric(df["starts"], errors='coerce').fillna(0).to_numpy()
	ends_arr = pd.to_numeric(df["ends"], errors='coerce').fillna(0).to_numpy()
	
	# Calculate width
	raw_width = ends_arr - starts_arr + 1
	
	# Force any value < 1 to become 1
	safe_width = np.clip(raw_width, a_min=1, a_max=None)
	
	gr = GenomicRanges(
		seqnames=df["seqnames"].tolist(),
		ranges=IRanges(
			start=starts_arr.astype(int).tolist(), 
			width=safe_width.astype(int).tolist()
			),
		strand=df["strand"].tolist() if "strand" in df.columns else ["*"] * len(df),
		mcols=BiocFrame(meta_df.to_dict("list"))
		)
	return gr

# Define gtf file paths
gtf_path = "data/gencode.v49.primary_assembly.annotation.gtf.gz"
parquet_path = "data/gencode.v49.primary_assembly.annotation.parquet"

# load GTF data
if os.path.exists(parquet_path):
	print(f"loading GTF data from Parquet file...")
	gtf_df = pd.read_parquet(parquet_path)
else:
	print("no cache found: parsing raw GTF (this will take a while)...")
	gtf_gr = genomicranges.read_gtf(gtf_path)
	gtf_df = gtf_gr.to_pandas()

	# Save it to Parquet for the next time
	print("Saving to Parquet for faster loading next time...")
	gtf_df.to_parquet(parquet_path, index=False)

# convert gene IDs to raw gene IDs without .1 or similar suffix
gtf_df['base_gene_id'] = gtf_df['gene_id'].astype(str).str.split('.').str[0]



################################
## 03 - make variant data     ##
## gr object and get pos info ##
################################

chunk_size = 5*10**6
first_chunk = True
total_variant_counts = pd.Series(dtype=int)

print(f"annotating variant positions (this may take a while)...")


with gzip.open(destination_filename, 'wt') as outfile:

	for input_df in pd.read_csv(input_filename, compression='gzip', sep='\t', chunksize=chunk_size, low_memory=False):
		
		# convert gene IDs to raw gene IDs without .1 or similar suffix
		input_df['base_gene_id'] = input_df['gene_id'].astype(str).str.split('.').str[0]

		# filter gtf genes only by target gene IDs
		target_gene_ids = input_df["base_gene_id"].tolist()
		mask = gtf_df["base_gene_id"].isin(target_gene_ids)
		target_gtf = gtf_df[mask]

		# Extract the full gene boundaries and exonic regions then convert to gr objects
		genes_df = target_gtf[target_gtf["feature"] == "gene"]
		exons_df = target_gtf[target_gtf["feature"] == "exon"]
		genes_gr = make_gr_from_df(genes_df)
		exons_gr = make_gr_from_df(exons_df)

		# make gr object from positive genes
		input_df = input_df.rename(columns={
			"chromosome": "seqnames", 
			"position": "starts"
			})
		input_df["starts"] = pd.to_numeric(input_df["starts"], errors='coerce')
		input_df["seqnames"] = input_df["seqnames"].astype(str)
		input_df["seqnames"] = input_df["seqnames"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
		input_df["ends"] = input_df["starts"] + input_df["ref"].str.len() - 1
		variants_gr = make_gr_from_df(input_df)

		# Calculate overlaps for both genes and exons
		exon_counts = exons_gr.count_overlaps(variants_gr)
		gene_counts = genes_gr.count_overlaps(variants_gr)

		# Loop through variant indeces and count exonic, intronic overlap
		classifications = []
		for e_count, g_count in zip(exon_counts, gene_counts):
			if e_count > 0:
				classifications.append("exonic")
			elif g_count > 0:
				classifications.append("intronic")
			else:
				classifications.append("intergenic")

		# Store resulting annotations in own column
		input_df["variant_location"] = classifications

		# revert 'seqnames' column back to 'chromosome'
		input_df = input_df.rename(columns={"seqnames": "chromosome"})
		input_df = input_df.drop(columns=["base_gene_id"], errors='ignore')

		# move 'ends' column behind 'starts'
		cols = list(input_df.columns)
		cols.remove("ends")
		starts_idx = cols.index("starts")
		cols.insert(starts_idx + 1, "ends")
		cols.remove("variant_location")
		starts_idx = cols.index("type")
		cols.insert(starts_idx + 1, "variant_location")
		input_df = input_df[cols]

		# Keep track of statistics across all chunks
		chunk_counts = input_df["variant_location"].value_counts()
		total_variant_counts = total_variant_counts.add(chunk_counts, fill_value=0)

		# save output to destination
		if first_chunk and not input_df.empty:
			input_df.to_csv(outfile, sep='\t', header=True, index=False)
			first_chunk = False
		elif not input_df.empty:
			input_df.to_csv(outfile, sep='\t', header=False, index=False)


total_time = time.time() - start_time 
print(f"Finished!\n {total_variant_counts.astype(int)}\nVariants processed in {total_time/60:.2f} minutes!")
print('')















