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

Step 1.: Use eQTL catalogue to identify variants

Each eQTL study in the catalog provides files on gene expression (ge) and also
on isoform usage (tx). The data is stored in .tsv.gz files and we want to extract
the following columns:
- "variant" 	variant id
- "gene_id"	id of the associated gene
- "pvalue"	significance of the association


