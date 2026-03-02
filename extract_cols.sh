#!/bin/bash


# This script extracts the following columns from a .tsv.gz file for storing eQTL data:
#- variant
#- gene_id
#- pvalue
#- molecular_trait_id
#- median_tpm
#and saves the column values in a user-defined output directory.


# get inputs
source_filename=$1
destination=$2

if [[ $source_filename == "" ]]; then
	echo "Missing source file name argument"
	exit 1
fi

if [[ $destination == "" ]]; then
	echo "Missing destination directory argument"
	exit 1
fi


# find column indeces
header=$(zcat "$source_filename" | head -1 | tr '\t' '\n' | nl)
col1=$(echo "$header" | grep 'variant' | cut -f1 | xargs)
col2=$(echo "$header" | grep 'gene_id' | cut -f1 | xargs)
col3=$(echo "$header" | grep 'pvalue' | cut -f1 | xargs)
col4=$(echo "$header" | grep 'molecular_trait_id' | cut -f1 | xargs)
col5=$(echo "$header" | grep 'median_tpm' | cut -f1 | xargs)

# Check if all found
if [[ -z "$col1" || -z "$col2" || -z "$col3" || -z "$col4" || -z "$col5" ]]; then
	echo "ERROR: Missing columns! Found: $col1,$col2,$col3,$col4,$col5"
	exit 1
fi

echo "Extracting cols $col1,$col2,$col3,$col4,$col5 to $destination"

# extract columns
zcat "$source_filename" | cut -f"$col1,$col2,$col3,$col4,$col5" | gzip > "$destination"
