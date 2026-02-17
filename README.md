# Special Course Project on Prediction of the Effects of eQTL Variants on Gene Expression or Isoform Usage

This is the repository for my special course project on prediction of 
differential gene expression or differential isoform usage based on eQTL data.

================================

The pipeline uses data from eQTL studies and the variant sequences as training
data for two ML classifiers. One to predict isoform usage, one to predict gene
expression.

The key steps are the following:
1. Identify variants significantly associated with a change in expression or isoform usage and negative controls.
2. Extract variant DNA sequences
3. Encode DNA sequences +200 NT window using a DNA language model
4. Build classifiers trained on embedded sequences

In the following, these steps are broken down a bit further and associated with 
the respective scripts in the repo.

================================

#Pipeline (run on node08):

### General pipeline
```shell
cd "/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction"

./extract_cols.bsh <eQTL_file.tsv.gz> <output_destination.tsv.gz>

python3 adjust_pv.py <extracted_cols.tsv.gz> <output_destination.tsv.gz> <alpha>

```

### Run example (Alasoo_2018 expression file)
```shell
cd "/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction"

./extract_cols.bsh /home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz results/Alasoo_2018_ge_extracted.tsv.gz

python3 adjust_pv.py results/Alasoo_2018_ge_extracted.tsv.gz results/Alasoo_2018_ge_adj.tsv.gz 0.05

```
