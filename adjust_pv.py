#!/bin/env/python3


##################################################
## this script bonferroni-corrects the pvalues 	##
##   in a .tsv.gz file containing eQTL data.	##
##################################################


# load libraries
import sys
import gzip
import pandas as pd
import numpy as np


# define functions
def adjust_pv(pv, m):
	p_adj = pv*m
	return min(p_adj, 1.0)


# store user input 
filename = sys.argv[1]		# name of .gz zipped file
pv_cutoff = sys.argv[2]		# user defined pv_cutoff


# input handling
if pv_cutoff == '' or filename == '':
	print("Error: missing input argument pvalue cutoff or file name!")
	sys.exit(1)
if not filename.endswith('.tsv.gz'):
	print("Error: input file should be .tsv.gz format")
	sys.exit(1)
try: 
	pv_cutoff = float(pv_cutoff)
except ValueError:
	print("Error: cutoff input should be float, e.g. 0.05")
	sys.exit(1)


# get m
m = 0
with gzip.open(filename, 'rt') as eQTL_file:
	for line in eQTL_file:
		m += 1
print(f'\nworking with m = {m}...')
print(f'adjusted significance level is alpha = {pv_cutoff/m}\n')


# load tsv
