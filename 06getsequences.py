#!/home/ctools/opt/anaconda-2025-12-2/bin/python3

import os
import sys
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
	print(f"python3 {sys.argv[0]} <positive_variant_filename> <negative_variant_filename> <reference_filename> <destination>")
	sys.exit(1)

# get user input for filenames
if sys.argv[1] == 'help':
	usage()

# handle filenames
try:
	pos_variant_filename = sys.argv[1]
	neg_variant_filename = sys.argv[2]
	reference_filename = sys.argv[3]
	destination_filename = sys.argv[4]
except IndexError:
	print(f"I am missing a filename.")
	usage()
if not pos_variant_filename.endswith('.tsv.gz') or not neg_variant_filename.endswith('.tsv.gz') or not destination_filename.endswith('.tsv.gz'):
	print(f"Make sure that variant and destination are .tsv.gz format!")
	print(f"\tpositives_filename = {pos_variant_filename}")
	print(f"\tpositives_filename = {neg_variant_filename}")
	print(f"\tdestination_filename = {destination_filename}")
	usage()
if not reference_filename.endswith('.fa.bgz'):
	print(f"Make sure that reference file is .fa.bgz format!")
	print(f"\tGiven: {reference_filename}")
	print(f"\tto convert .fa.gz to .fa.bgz you can run for example:")
	print("\t\tgunzip -c data/GRCh38.primary_assembly.genome.fa.gz | bgzip > data/GRCh38.primary_assembly.genome.fa.bgz")
	usage()
if not os.path.exists(pos_variant_filename): 
	print(f"{pos_variant_filename} was not found at specified location.")
	sys.exit(1)
elif not os.path.exists(neg_variant_filename): 
	print(f"{neg_variant_filename} was not found at specified location.")
	sys.exit(1)
elif not os.path.exists(reference_filename): 
	print(f"{reference_filename} was not found at specified location.")
	sys.exit(1)




##############################
##   02 - Get sequence str  ##
##############################


# function to extract sequence
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
	start_idx = row['starts'] - 1 - window
	stop_idx = row['ends'] + window

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

# get variant file
print(f"Reading variant files...")
posvar_df = pd.read_csv(pos_variant_filename, compression='gzip', sep='\t', low_memory=False)
negvar_df = pd.read_csv(neg_variant_filename, compression='gzip', sep='\t', low_memory=False)

# extract sequences and store in a new column
print(f"extracting sequences...")
posvar_df['variant_window'] = posvar_df.apply(lambda row: get_seq(row, gene_ref=all_genes, window=100), axis=1)
negvar_df['variant_window'] = negvar_df.apply(lambda row: get_seq(row, gene_ref=all_genes, window=100), axis=1)

# store as .tsv.gz
print(f"storing results...")
combined_df = pd.concat([posvar_df, negvar_df], ignore_index=True)
combined_df.to_csv(destination_filename, sep='\t', index=False, compression='gzip')

total_time = time.time() - start_time
print(f'finished after {total_time/60:.2f} minutes!') 
print(f"results can be found at {destination_filename}\n")






