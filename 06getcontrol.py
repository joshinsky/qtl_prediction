#!/home/ctools/opt/anaconda-2025-12-2/bin/python3

import os
import sys
import pandas as pd
import numpy as np
import time
start_time = time.time()


##############################
##   01 - Get user input    ##
## (filenames and criteria) ##
##############################

# function explaining the program usage
def usage():
	print("program usage:")
	print(f"python3 {sys.argv[0]} <positives_filename> <non_sig_filename> <destination> [criteria]")
	print(f"criteria could be: gene, location, variant")
	print(f"the order of criteria determines their importance\n")
	sys.exit(1)

# get user input for filenames
if sys.argv[1] == 'help':
	usage()

# handle filenames
try:
	positives_filename = sys.argv[1]
	non_sig_filename = sys.argv[2]
	destination_filename = sys.argv[3]
except IndexError:
	print(f"I am missing a filename.")
	usage()

if not positives_filename.endswith('.tsv.gz') or not non_sig_filename.endswith('.tsv.gz') or not destination_filename.endswith('.tsv.gz'):
	print(f"I am missing a filename. These were given:")
	print(f"\tpositives_filename = {positives_filename}")
	print(f"\tnon_sig_filename = {non_sig_filename}")
	print(f"\tdestination_filename = {destination_filename}")
	print(f"Make sure that input is .tsv.gz format!")
	usage()
if not os.path.exists(positives_filename): 
	print(f"{positives_filename} was not found at specified location.")
	sys.exit(1)
elif not os.path.exists(non_sig_filename): 
	print(f"{non_sig_filename} was not found at specified location.")
	sys.exit(1)
else:
	print(f"Reading filters.")


# get user input for filtering criteria
criteria_args = sys.argv[4:]
if len(criteria_args) == 0:
	print(f"No criteria args were input.")
	usage()


################################
##     02 - Get criteria      ##
## (extract+rank target cols) ##
################################

# list available criteria
criterion_to_col = {"gene":		"gene_id",
					"location":	"variant_location",
					"variant":	"type"}

# get criteria to use in a ranked list
criteria_used = []
for c in criteria_args:
	if c in criterion_to_col:
		criteria_used.append(criterion_to_col[c])

# handle unusable/invalid criteria
if len(criteria_used) == 0:
	print(f"input criteria were not valid.")
	usage()
if len(criteria_used) != len(criteria_args):
	print(f"warning! Could not find all desired criteria.")
	print(f"input: {criteria_args}")
	print(f"using: {criteria_used}")
	print(f"to interrupt, press ctrl+C")
else:
	print(f"using: {criteria_used}")


#########################################
##          03 - Get negatives         ##
## (find most similar non sig entries) ##
#########################################

# load data
print(f"loading data...")
pos_df = pd.read_csv(positives_filename, compression='gzip', sep='\t', low_memory=False)
nsig_df = pd.read_csv(non_sig_filename, compression='gzip', sep='\t', low_memory=False)

print(f"finding negative control candidates...")
neg_df = pd.merge(pos_df, nsig_df, how='inner', on=criteria_used)

if len(neg_df.index) < len(pos_df.index):
	print("could not find negative control for all positive controls!")
	print(f"found {len(neg_df.index)} negatives for {len(pos_df.index)} positives.")

# keep only the first gene in case of multiple perfect matches
print(f"selecting negative controls...")
neg_df = neg_df.drop_duplicates(subset=['gene_id'])

print(f"storing results...")
# control_df = pd.DataFrame(control_rows)
neg_df.to_csv(destination_filename, sep='\t', index=False, compression='gzip')

total_time = time.time() - start_time
print(f'finished after {total_time/60:.2f} minutes!') 
print(f"results can be found at {destination_filename}\n")

