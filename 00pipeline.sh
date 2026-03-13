#!/bin/bash

# stop pipeline on error
set -e
set -o pipefail

# get user input
location=$1
test=$2

################
## get set up ##
################

# setting up environment
if [ ! -d "data" ]; then
	echo "making /data directory."
	mkdir data
fi

if [ ! -d "temp" ]; then
	echo "making /temp directory."
	mkdir temp
fi

if [ ! -d "results" ]; then
	echo "making /results directory."
	mkdir results
fi

# get test file
if [[ $location == "cluster" ]]; then
	cd /home/projects2/kvs_students/2026/jl_qtl_prediction
	source /home/ctools/opt/anaconda-2025-12-2/bin/activate josh_env
	cd /home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction
	test_file=/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz
elif [[ $location == "josh" ]]; then
	cd "/Users/Study/Library/CloudStorage/OneDrive-Personal/002 Education/002 Uni/004 DTU/004 Semester/26_02_S2/specou/repo"
	test_file=data/test_Alasoo2018.tsv.gz
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
	echo "running test for $test_file on first 5 Mio. rows"
	echo
	gunzip -c $test_file | head -5000001 | gzip > temp/test_subset.tsv.gz
	startfile=temp/test_subset.tsv.gz
else
	echo
	echo
	echo "running pipeline for $test_file"
	echo
	startfile=$test_file
fi


##############################
## Bonferoni correct pvalue ##
##############################

# adjust pvalue and split in sig and non sig
echo
echo "run 01adjpv.py"
if [[ $location == "cluster" ]]; then
	./01adjpv.py $startfile temp/test_padj.tsv.gz 0.05 0.9
elif [[ $location == "josh" ]]; then
	python3 01adjpv.py $startfile temp/test_padj.tsv.gz 0.05 0.9
fi


######################
## get positive set ##
######################

# split on sig 
echo
echo "run 02splitsig.py"
if [[ $location == "cluster" ]]; then
	./02splitsig.py temp/test_padj.tsv.gz temp/test_sig.tsv.gz temp/test_nonsig.tsv.gz 
elif [[ $location == "josh" ]]; then
	python3 02splitsig.py temp/test_padj.tsv.gz temp/test_sig.tsv.gz temp/test_nonsig.tsv.gz 
fi

# get most sig variants per gene
echo
echo "run 03topsig.py"
if [[ $location == "cluster" ]]; then
	./03topsig.py temp/test_sig.tsv.gz temp/test_most_sig.tsv.gz
elif [[ $location == "josh" ]]; then
	python3 03topsig.py temp/test_sig.tsv.gz temp/test_most_sig.tsv.gz
fi


##########################
## Add positional info  ##
##########################

# retrieve data
if [ ! -f "data/gencode.v49.primary_assembly.annotation.gtf.gz" ]; then
	echo
	echo "retrieving annotation file for positional info"
	cd data
	wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
	cd ..
fi

# add positional info for sigs
echo
echo "run 04getposition.py on significant variants"
if [[ $location == "cluster" ]]; then
	./04getposition.py temp/test_most_sig.tsv.gz temp/test_positives.tsv.gz
elif [[ $location == "josh" ]]; then
	python3 04getposition.py temp/test_most_sig.tsv.gz temp/test_positives.tsv.gz
fi

# add positional info for non-sigs
echo
echo "run 04getposition.py on non-significant variants"
if [[ $location == "cluster" ]]; then
	./04getposition.py temp/test_nonsig.tsv.gz temp/test_nonsig_annotated.tsv.gz
elif [[ $location == "josh" ]]; then
	python3 04getposition.py temp/test_nonsig.tsv.gz temp/test_nonsig_annotated.tsv.gz
fi


###########################
## Get negative controls ##
###########################

# get negative controls
echo
echo "run 05getnegatives.py"
if [[ $location == "cluster" ]]; then
	./05getnegatives.py temp/test_positives.tsv.gz temp/test_nonsig_annotated.tsv.gz temp/test_negatives.tsv.gz gene location variant
elif [[ $location == "josh" ]]; then
	python3 05getnegatives.py temp/test_positives.tsv.gz temp/test_nonsig_annotated.tsv.gz temp/test_negatives.tsv.gz gene location variant
fi


#######################
## add sequence data ##
#######################

# retrieve data
if [ ! -f "data/GRCh38.primary_assembly.genome.fa.bgz.fai" ]; then
	if [ ! -f "data/GRCh38.primary_assembly.genome.fa.gz" ]; then
		cd data
		echo
		echo "retrieving reference genome"
		wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
		cd ..
	fi
	echo
	echo "converting reference genome to bgzip for pyfaidx to use."
	gunzip -c data/GRCh38.primary_assembly.genome.fa.gz | bgzip > data/GRCh38.primary_assembly.genome.fa.bgz
fi

echo
echo "run 06getsequences.py"
if [[ $location == "cluster" ]]; then
	./06getsequences.py temp/test_positives.tsv.gz temp/test_negatives.tsv.gz data/GRCh38.primary_assembly.genome.fa.bgz results/test_seqs.tsv.gz
elif [[ $location == "josh" ]]; then
	python3 06getsequences.py temp/test_positives.tsv.gz temp/test_negatives.tsv.gz data/GRCh38.primary_assembly.genome.fa.bgz results/test_seqs.tsv.gz
fi

echo
echo "############################################"
echo "finished for $test_file"
echo "############################################"
echo
echo