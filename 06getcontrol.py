#!/home/ctools/opt/anaconda-2025-12-2/bin/python3

import os
import sys
import pandas as pd
import numpy as np


##############################
##   01 - Get user input    ##
## (filenames and criteria) ##
##############################

# function explaining the program usage
def usage():
	print("\nprogram usage:")
	print(f"python3 {sys.argv[0]} <positives_filename> <non_sig_filename> <destination> [criteria]")
	print(f"\ncriteria could be: gene, position, variant")
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
	print(f"\nI am missing a filename.")
	usage()

if not positives_filename.endswith('.tsv.gz') or not non_sig_filename.endswith('.tsv.gz') or not destination_filename.endswith('.tsv.gz'):
	print(f"\nI am missing a filename. These were given:")
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
criterion_to_col = {
    "gene":     "gene_id",
    "position": "rel_pos",
    "variant":     "type"
}

# get criteria to use in a ranked list
criteria_used = []
for c in criteria_args:
    if c in criterion_to_col:
        criteria_used.append(criterion_to_col[c])

# handle unusable/invalid criteria
if len(criteria_used) == 0:
	print(f"Input criteria were not valid.")
	usage()
if len(criteria_used) != len(criteria_args):
	print(f"Warning! Could not find all desired criteria.")
	print(f"Input: {criteria_args}")
	print(f"Using: {criteria_used}")
	print(f"To interrupt, press ctrl+C")
else:
	print(f"Using: {criteria_used}\n")


# #########################################
# ##          03 - Get negatives         ##
# ## (find most similar non sig entries) ##
# #########################################

# # load data
# print(f"loading data...\n")
# pos_df = pd.read_csv(positives_filename, compression='gzip', sep='\t')
# nsig_df = pd.read_csv(non_sig_filename, compression='gzip', sep='\t')

# # initialise loop
# control_rows = []  
# criteria_cols = [criterion_to_col[crit] for crit in criteria_used]

# # for each significant entry find the best (most similar) negative control
# print(f"finding negative controls for {len(pos_df.index)} significant entries...\n")
# for pos_idx in range(len(pos_df)):
#     pos_row = pos_df.iloc[pos_idx]
    
#     best_match_vector = [-1] * len(criteria_used)
#     best_candidates = []
    
#     for neg_idx in range(len(nsig_df)):
#         neg_row = nsig_df.iloc[neg_idx]

#         # get match between columns as a binary vector
#         match_vector = [1 if pos_row[col] == neg_row[col] else 0 for col in criteria_cols]
        
#         # test if current match is better than current best match
#         if match_vector > best_match_vector:
#             best_match_vector = match_vector[:]
#             best_candidates = [neg_row]
#         elif match_vector == best_match_vector:
#             best_candidates.append(neg_row)
    
#     # Pick first best candidate
#     control = best_candidates[0]
#     control_rows.append(control)
#     print(f"found match for {pos_row['gene_id']} --> matching vector = {best_match_vector}:")

# print(f"Storing results")
# control_df = pd.DataFrame(control_rows)
# control_df.to_csv(destination_filename, sep='\t', index=False, compression='gzip')

# print(f"Finished! results can be found at {destination_filename}")
# print(f"Confirm using:\ngunzip -c {destination_filename} | head")



#########################################
##          03 - Get negatives         ##
## (find most similar non sig entries) ##
#########################################


# load data
print(f"loading data...\n")
pos_df = pd.read_csv(positives_filename, compression='gzip', sep='\t')
nsig_df = pd.read_csv(non_sig_filename, compression='gzip', sep='\t')

negatives_df = pd.merge(pos_df, nsig_df, how='inner', on=criteria_used)







































