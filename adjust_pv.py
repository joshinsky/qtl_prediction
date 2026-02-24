#!/bin/env/python3


##################################################
## this script bonferroni-corrects the pvalues 	##
##   in a .tsv.gz file containing eQTL data.	##
##################################################


# load libraries
import sys
import gzip
import math
import pandas as pd
import numpy as np
import time
start_time = time.time()


# store user input 
filename = sys.argv[1]          # name of .gz zipped file
destination = sys.argv[2]	# name of .gz zipped file for results
pv_cutoff = sys.argv[3]         # user defined pv_cutoff
pv_cutoff_non_sig = sys.argv[4]	# user defined pv_cutoff for non-significant

# input handling
if pv_cutoff == '' or pv_cutoff_non_sig == '' or filename == '' or destination == '':
	print("Error: missing input argument!")
	print("Usage:\npython3 adjust_pv.py <filename> <destination> <pv_cutoff> <pv_cutoff_non_sig>")
	sys.exit(1)
if not filename.endswith('.tsv.gz') or not destination.endswith('.tsv.gz'):
	print("Error: input file and destination should be in .tsv.gz format")
	sys.exit(1)
try: 
	pv_cutoff = float(pv_cutoff)
	pv_cutoff_non_sig = float(pv_cutoff_non_sig)
except ValueError:
	print("Error: pv_cutoff input should be float, e.g. 0.05")
	sys.exit(1)


# get m = number of lines
print(f'\ncalculating adjusted cutoff...')
m = -1				# subtract header line
with gzip.open(filename, 'rt') as eQTL_file:
        for line in eQTL_file:
                m += 1
print(f'\nworking with 	m = {m}...')
print(f'        alpha_adj = {pv_cutoff/m}...')
print(f'    non_sig_alpha = {pv_cutoff_non_sig}')


# prepare pvalue adjustment, define necessary variables
print(f'\nadjusting p-values...')
total_chunks = math.ceil(m/(5*10**6))
first_chunk = True
chunk_num = 0
sig_tot = 0

# get significant and non-significant entries as bool (1|0)
for chunk in pd.read_csv(filename, compression='gzip', sep='\t', usecols=['variant', 'gene_id', 'pvalue'], dtype={'variant':'string', 'gene_id':'string', 'pvalue':'float32'}, chunksize=5*10**6):
	p_value = chunk['pvalue']
	p_adj = np.minimum(p_value*m, 1.0)

	# store p_adj, sig and non_sig in .tsv.gz
	chunk['p_adj'] = p_adj
	chunk['significant'] = (p_adj <= pv_cutoff).astype(int)
	chunk['non_significant'] = (p_value >= pv_cutoff_non_sig).astype(int)
	chunk.to_csv(destination, mode='a', header=first_chunk, compression='gzip', index=False)
	
	first_chunk = False
	chunk_num +=1

	# get stats for continuous output while running
	sig = chunk['significant'].sum()
	non_sig = chunk['non_significant'].sum()
	sig_tot += sig
	non_sig_tot += non_sig
	print(f"Chunk {chunk_num}/{total_chunks}, \nsignificant: {sig}, total: {sig_tot} \nnon-significant: {non_sig}, total: {non_sig_tot}")


# final message
total_time = time.time() - start_time
print(f'analysis finished after {total_time/60:.2f} minutes!') 
print(f'confirm success using:')
print(f'zcat {destination} | head -25')
