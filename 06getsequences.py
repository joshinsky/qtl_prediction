#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import os
import sys
import gzip
import pandas as pd
import numpy as np
import time
from pyfaidx import Fasta
start_time = time.time()


##############################
##   01 - Get user input    ##
## (filenames and criteria) ##
##############################

# function explaining the program usage
def usage():
	print("program usage:")
	print(f"python3 {sys.argv[0]} <variant_filename> <reference_filename> <destination>")
	sys.exit(1)

# get user input for filenames
if sys.argv[1] == 'help':
	usage()

# handle filenames
try:
	variant_filename = sys.argv[1]
	reference_filename = sys.argv[2]
	destination_filename = sys.argv[3]
except IndexError:
	print(f"I am missing a filename.")
	usage()
if not variant_filename.endswith('.tsv.gz') or not destination_filename.endswith('.tsv.gz'):
	print(f"Make sure that variant and destination are .tsv.gz format!")
	print(f"\tvariant_filename = {variant_filename}")
	print(f"\tdestination_filename = {destination_filename}")
	usage()
if not reference_filename.endswith('.fa.bgz'):
	print(f"Make sure that reference file is .fa.bgz format!")
	print(f"\tGiven: {reference_filename}")
	print(f"\tto convert .fa.gz to .fa.bgz you can run for example:")
	print("\t\tgunzip -c data/GRCh38.primary_assembly.genome.fa.gz | bgzip > data/GRCh38.primary_assembly.genome.fa.bgz")
	usage()
if not os.path.exists(variant_filename): 
	print(f"{variant_filename} was not found at specified location.")
	sys.exit(1)
elif not os.path.exists(reference_filename): 
	print(f"{reference_filename} was not found at specified location.")
	sys.exit(1)




##############################
##   02 - Get sequence str  ##
##############################

# helper function to extract sequence window
def get_seq(row, gene_ref, window=0):
	"""	This function extracts the DNA sequence of a section on the chromosome using 
		start and stop indeces.
	"""

	chromosome = row['chromosome']
	if chromosome not in gene_ref.keys():
		print(f"error: reference does not contain chromosome '{chromosome}'.")
		sys.exit(1)
	if not isinstance(window, int):
		print(f"error: nucleotide window should be given as integer not as {type(window)}.")
		sys.exit(1)

	# get start and end indeces, correct for 1-based indexing in pyfaidx
	start = int(row['starts'])
	stop = int(row['ends'])
	chrom_len = len(gene_ref[chromosome])
	start_idx = max(0, start - 1 - window)
	stop_idx = min(chrom_len, stop + window)

	# extract sequence as plain string
	sequence = gene_ref[chromosome][start_idx:stop_idx].seq

	if sequence == "":
		print(f"error: no sequence found for {row['gene_id']}.")
		return pd.NA
	else:
		return sequence

# get reference sequences
print(f"Reading genome reference...")
if not os.path.exists(reference_filename+'.fai'):
	print(f"storing indexed file {reference_filename+'.fai'} for later fast lookups.")
all_genes = Fasta(reference_filename)

# open variant file in chunks, then annotate with sequence
print(f"extracting sequence for each variant...")
chunk_size = 10**6
first_chunk = True
with gzip.open(destination_filename, 'wt') as outfile:
	
	for df_chunk in pd.read_csv(variant_filename, compression='gzip', sep='\t', chunksize=chunk_size, low_memory=False):
		
		# extract sequences and store in a new column
		df_chunk['variant_window'] = df_chunk.apply(lambda row: get_seq(row, gene_ref=all_genes, window=100), axis=1)

		# write to .tsv.gz
		if first_chunk and not df_chunk.empty:
			df_chunk.to_csv(outfile, sep='\t', header=True, index=False)
			first_chunk = False
		elif not df_chunk.empty:
			df_chunk.to_csv(outfile, sep='\t', header=False, index=False)

# final message
total_time = time.time() - start_time
print(f'finished after {total_time/60:.2f} minutes!') 
print(f"results can be found at {destination_filename}\n")






