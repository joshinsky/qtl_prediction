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
	print(f"\tpositives_filename = {variant_filename}")
	print(f"\tdestination_filename = {destination_filename}")
	usage()
if not reference_filename.endswith('.fa.gz'):
	print(f"Make sure that reference file is .fa.gz format!")
	print(f"\tGiven: {reference_filename}")
	usage()

if not os.path.exists(variant_filename): 
	print(f"{variant_filename} was not found at specified location.")
	sys.exit(1)
elif not os.path.exists(reference_filename): 
	print(f"{reference_filename} was not found at specified location.")
	sys.exit(1)
else:
	print(f"Reading genome reference.")




##############################
##   02 - Read genome ref   ##
##                          ##
##############################













