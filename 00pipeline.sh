#!/bin/bash

location=$1
test=$2

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

# subset to 5,000,000 rows + header
if [[ $test == "test" ]]; then
	echo
	echo
	echo "running test for $raw_file on first 5 Mio. rows"
	echo
	gunzip -c $raw_file | head -5000001 | gzip > temp/test_subset.tsv.gz
	startfile=temp/test_subset.tsv.gz
else
	echo
	echo
	echo "running pipeline for $raw_file"
	echo
	startfile=$raw_file
fi

# subset columns
# ./01getcols.sh temp/test_subset.tsv.gz temp/test_subset_cols.tsv.gz

# adjust pvalue and split in sig and non sig
echo
echo "run 02adjpv.py"
if [[ $location == "cluster" ]]; then
	./02adjpv.py $startfile temp/test_padj.tsv.gz 0.05 0.9
elif [[ $location == "josh" ]]; then
	python3 02adjpv.py $startfile temp/test_padj.tsv.gz 0.05 0.9
fi

# split on sig 
echo
echo "run 03splitsig.py"
if [[ $location == "cluster" ]]; then
	./03splitsig.py temp/test_padj.tsv.gz temp/test_sig.tsv.gz temp/test_nonsig.tsv.gz 
elif [[ $location == "josh" ]]; then
	python3 03splitsig.py temp/test_padj.tsv.gz temp/test_sig.tsv.gz temp/test_nonsig.tsv.gz 
fi

# get most sig variants per gene
echo
echo "run 04topsig.py"
if [[ $location == "cluster" ]]; then
	./04topsig.py temp/test_sig.tsv.gz temp/test_most_sig.tsv.gz
elif [[ $location == "josh" ]]; then
	python3 04topsig.py temp/test_sig.tsv.gz temp/test_most_sig.tsv.gz
fi

# add positional info for sigs
echo
echo "run 05getpos.py on significant variants"
if [[ $location == "cluster" ]]; then
	./05getpos.py temp/test_most_sig.tsv.gz temp/test_positives.tsv.gz
elif [[ $location == "josh" ]]; then
	python3 05getpos.py temp/test_most_sig.tsv.gz temp/test_positives.tsv.gz
fi

# add positional info for non-sigs
echo
echo "run 05getpos.py on non-significant variants"
if [[ $location == "cluster" ]]; then
	./05getpos.py temp/test_nonsig.tsv.gz temp/test_nonsig_annotated.tsv.gz
elif [[ $location == "josh" ]]; then
	python3 05getpos.py temp/test_nonsig.tsv.gz temp/test_nonsig_annotated.tsv.gz
fi

# get negative controls
echo
echo "run 06getcontrol.py"
if [[ $location == "cluster" ]]; then
	./06getcontrol.py temp/test_positives.tsv.gz temp/test_nonsig_annotated.tsv.gz temp/test_negatives.tsv.gz gene location variant
elif [[ $location == "josh" ]]; then
	python3 06getcontrol.py temp/test_positives.tsv.gz temp/test_nonsig_annotated.tsv.gz temp/test_negatives.tsv.gz gene location variant
fi


echo
echo "############################################"
echo "finished for $raw_file"
echo "############################################"
echo
echo