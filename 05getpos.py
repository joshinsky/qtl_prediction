#!/home/ctools/opt/anaconda-2025-12-2/bin/python3

import os
import sys
import pandas as pd
import numpy as np
import genomicranges
from genomicranges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
# import biostrings

##############################
##   01 - Get user input    ##
## (filenames)              ##
##############################

# function explaining the program usage
def usage():
	print("\nprogram usage:")
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
	print(f"\nI am missing a filename.")
	usage()

##############################
## 02 - Load gtf data and   ##
##  get gr objects          ##
##############################

# function for translating a data frame into a gr object
def make_gr_from_df(df):
	meta_df = df.drop(columns=["seqnames", "starts", "ends", "strand"], errors="ignore")
	gr = GenomicRanges(
		seqnames=df["seqnames"].tolist(),
		ranges=IRanges(
			start=df["starts"].tolist(), 
			width=(df["ends"] - df["starts"]).astype(int).tolist()
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
	print(f"loading GTF data from Parquet file...\n")
	gtf_df = pd.read_parquet(parquet_path)
else:
	print("No cache found. Parsing raw GTF (this will take a while)...\n")
	gtf_gr = genomicranges.read_gtf(gtf_path)
	gtf_df = gtf_gr.to_pandas()

	# Save it to Parquet for the next time
	print("Saving to Parquet for faster loading next time...\n")
	gtf_df.to_parquet(parquet_path, index=False)

# load gene variant data
print(f"loading variant data...\n")
input_df = pd.read_csv(input_filename, compression='gzip', sep='\t')

# filter gtf genes only by target gene IDs
print(f"extracting target genes...\n")

gtf_df['base_gene_id'] = gtf_df['gene_id'].astype(str).str.split('.').str[0]	# look only at raw gene IDs without .1 or similar suffix
input_df['base_gene_id'] = input_df['gene_id'].astype(str).str.split('.').str[0]
target_gene_ids = input_df["base_gene_id"].tolist()
mask = gtf_df["base_gene_id"].isin(target_gene_ids)
target_gtf = gtf_df[mask]

print(f"found {len(target_gtf)}/{len(input_df)} target genes.\n")

# Extract the full gene boundaries and exonic regions then convert to gr objects
print(f"retrieving positional indeces...\n")
genes_df = target_gtf[target_gtf["feature"] == "gene"]
exons_df = target_gtf[target_gtf["feature"] == "exon"]
genes_gr = make_gr_from_df(genes_df)
exons_gr = make_gr_from_df(exons_df)

# make gr object from positive genes
input_df = input_df.rename(columns={
	"chromosome": "seqnames", 
	"position": "starts"
	})
input_df["seqnames"] = input_df["seqnames"].astype(str)
input_df["seqnames"] = input_df["seqnames"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
input_df["ends"] = input_df["starts"] + input_df["ref"].str.len() - 1
variants_gr = make_gr_from_df(input_df)



##############################
## 03 - extract variant pos ##
## (intrnc, exonc, intergnc)##
##############################

# Calculate overlaps for both genes and exons
print(f"classifying positional overlaps for each variant...\n")
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

# Make column for resulting annotations
input_df["variant_location"] = classifications

# Save the final output
print(f"\nSaving results to {destination_filename}...\n")
input_df.to_csv(destination_filename, sep='\t', index=False, compression='gzip')

print("Summary:")
print(input_df["variant_location"].value_counts())
















