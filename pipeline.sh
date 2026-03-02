#!/bin/bash

mydir
cd repo/qtl_prediction

# subset to 5,000,000 rows (5%)
zcat /home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz | head -2000001 | gzip > temp/test_subset.tsv.gz

# subset columns
time bash extract_cols.sh temp/test_subset.tsv.gz temp/test_subset_cols.tsv.gz

# view results
# zcat temp/test_subset_cols.tsv.gz | head -15

# adjust pvalue and split in sig and non sig
python3 adjust_pv.py temp/test_subset_cols.tsv.gz temp/test_padj.tsv.gz 0.05 0.9

# view results
# zcat temp/test_padj.tsv.gz | head -25

# add sequence and positional data

# split on sig 
python3 split_sig.py temp/test_padj.tsv.gz temp/test_sig.tsv.gz temp/test_nonsig.tsv.gz 

zcat temp/test_sig.tsv.gz | cut -f5,6 | head -25
echo
zcat temp/test_nonsig.tsv.gz | cut -f5,6 | head -25