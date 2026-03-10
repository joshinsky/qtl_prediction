#!/bin/bash

location=$1

if [[ $location == "cluster" ]]; then
	cd /home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction
	raw_file=/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz
elif [[ $location == "josh" ]]; then
	cd "/Users/Study/Library/CloudStorage/OneDrive-Personal/002 Education/002 Uni/004 DTU/004 Semester/26_02_S2/specou/repo"
	raw_file=data/test_Alasoo2018
else
	echo "Please specify which machine you're working on:"
	echo "on healthtech cluster?-->	./pipeline.sh cluster "
	echo "on Joshua's laptop? 	-->	./pipeline.sh josh "
	exit 1
fi

# subset to 5,000,000 rows (5%) + header
gunzip -c $raw_file | head -5000001 | gzip > temp/test_subset.tsv.gz

# subset columns
./01getcols.sh temp/test_subset.tsv.gz temp/test_subset_cols.tsv.gz

# view results
# gunzip -c temp/test_subset_cols.tsv.gz | head -15

# adjust pvalue and split in sig and non sig
./02adjpv.py temp/test_subset_cols.tsv.gz temp/test_padj.tsv.gz 0.05 0.9

# view results
# gunzip -c temp/test_padj.tsv.gz | head -25

# add sequence and positional data

# split on sig 
./04splitsig.py temp/test_padj.tsv.gz temp/test_sig.tsv.gz temp/test_nonsig.tsv.gz 

# gunzip -c temp/test_sig.tsv.gz | cut -f5,6 | head -25
# echo
# gunzip -c temp/test_nonsig.tsv.gz | cut -f5,6 | head -25

# get most sig variants per gene
./05topsig.py temp/test_sig.tsv.gz temp/test_most_sig.tsv.gz

# get negative controls
./06getcontrol.py temp/test_most_sig.tsv.gz temp/test_nonsig.tsv.gz temp/test_negatives.tsv.gz gene variant