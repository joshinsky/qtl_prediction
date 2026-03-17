#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import sys
import gzip
import numpy as np
import pandas as pd
import time
start_time = time.time()

# store user input
source_filename = sys.argv[1]
sig_destination = sys.argv[2]
nsig_destination = sys.argv[3]

# input handling
if source_filename == '' or sig_destination == '' or nsig_destination == '':
	print("Error: missing input argument!")
	print("Usage:\npython3 split_sig.py <filename> <sig_destination> <non_sig_destination>")
	sys.exit(1)
if not source_filename.endswith('.tsv.gz') or not sig_destination.endswith('.tsv.gz') or not nsig_destination.endswith('.tsv.gz'):
	print("Error: input file and destination should be in .tsv.gz format")
	sys.exit(1)

# split file on sig and non_sig in chunks for memory
print(f"splitting {source_filename} on significance.")
chunk_size = 5*10**6
first_chunk = True
for chunk in pd.read_csv(source_filename, compression='gzip', sep='\t', chunksize=chunk_size, low_memory=False):
	# define split masks
	sig_mask = (chunk['significant'] == 1)
	nsig_mask = (chunk['non_significant'] == 1)

	# split and store
	sig_chunk = chunk[sig_mask]
	nsig_chunk = chunk[nsig_mask]
	if first_chunk:
		if not sig_chunk.empty:
			sig_chunk.to_csv(sig_destination, sep='\t', mode='w', header=True, compression='gzip', index=False)
		if not nsig_chunk.empty:
			nsig_chunk.to_csv(nsig_destination, sep='\t', mode='w', header=True, compression='gzip', index=False)
		first_chunk = False
	else:
		sig_chunk.to_csv(sig_destination, sep='\t', mode='a', header=False, compression='gzip', index=False)
		nsig_chunk.to_csv(nsig_destination, sep='\t', mode='a', header=False, compression='gzip', index=False)


# final message
total_time = time.time() - start_time
print(f'finished after {total_time/60:.2f} minutes!') 
print(f'significant set stored at: {sig_destination}')
print(f'non-significant set stored at: {nsig_destination}\n')
