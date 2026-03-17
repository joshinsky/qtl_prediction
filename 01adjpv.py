#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3


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
import subprocess
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
print(f'calculating adjusted pvalue cutoff for {filename}...')
result = subprocess.run(f"gunzip -c {filename} | wc -l", shell=True, capture_output=True, text=True)
m = int(result.stdout.strip()) - 1 	# subtract header line
print(f'working with 	m = {m}')
print(f'        alpha_adj = {pv_cutoff/m}')
print(f'    non_sig_alpha = {pv_cutoff_non_sig}')


# prepare pvalue adjustment, define necessary variables
if pv_cutoff >= pv_cutoff_non_sig:
	print(f"WARNING! You specified a significance cutoff larger than the non-significance cutoff:")
	print(f"alpha = {pv_cutoff}")
	print(f"non_sig_alpha = {pv_cutoff_non_sig}")
	print(f"perhaps you got the inputs mixed up? To abort press ctrl+C")
total_chunks = math.ceil(m/(5*10**6))
first_chunk = True
chunk_num = 0
sig_tot = 0
non_sig_tot = 0

# get significant and non-significant entries as bool (1|0)
print(f'Bonferoni correcting p-values in {total_chunks} chunks...')
with gzip.open(destination, 'wt') as outfile:
	for chunk in pd.read_csv(filename, compression='gzip', sep='\t', dtype={'variant':'string', 'gene_id':'string', 'pvalue':'float32'}, chunksize=5*10**6, low_memory=False):
		
		# subset to only necessary columns
		chunk.drop(columns=["ma_samples", "maf", "beta", "se", "ac", "an", "r2", "rsid"], inplace=True, errors='ignore')

		# adjust pvalue
		p_value = chunk['pvalue']
		p_adj = np.minimum(p_value*m, 1.0)
	
		# remove non-significant entries >alpha_non_sig
		mask = (p_adj <= pv_cutoff) | (p_value >= pv_cutoff_non_sig)
		filtered_chunk = chunk.loc[mask].copy()
	
		# store p_adj, sig and non_sig as new columns
		filtered_chunk['p_adj'] = p_adj[mask]
		filtered_chunk['significant'] = (filtered_chunk['p_adj'] <= pv_cutoff).astype(int)
		filtered_chunk['non_significant'] = (filtered_chunk['pvalue'] >= pv_cutoff_non_sig).astype(int)
		
		# write to .tsv.gz
		if first_chunk and not filtered_chunk.empty:
			filtered_chunk.to_csv(outfile, sep='\t', header=True, index=False)
			first_chunk = False
		elif not filtered_chunk.empty:
			filtered_chunk.to_csv(outfile, sep='\t', header=False, index=False)
		chunk_num +=1
	
		# get stats for continuous output while running
		sig = filtered_chunk['significant'].sum()
		non_sig = filtered_chunk['non_significant'].sum()
		sig_tot += sig
		non_sig_tot += non_sig
		print(f"chunk {chunk_num}/{total_chunks}, significant: {sig}, non-significant: {non_sig}")


# final message
total_time = time.time() - start_time
print(f'finished after {total_time/60:.2f} minutes!')
print(f'extracted {sig_tot} sig and {non_sig_tot} nonsig entries...')
print(f'stored at {destination}\n')

