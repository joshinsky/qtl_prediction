# Special Course Project on Prediction of eQTL Variant Effects on Gene Expression and Isoform Usage

This is the repository for my special course project on prediction of 
differential gene expression or differential isoform usage based on eQTL data.

================================

The pipeline uses data from eQTL studies and the variant sequences as training
data for two ML classifiers. One to predict isoform usage, one to predict gene
expression.

The key steps are the following:
1. Identify variants significantly associated with a change in expression or isoform usage
2. Find suitable negative controls (one set of closely related nsig variants and one set of random negatives)
3. Extract variant position on the chromosome and DNA sequences +-100 nucleotide window
4. Encode DNA sequences using a DNA language model
5. Build classifiers trained on embedded sequences for effect prediction

================================

Current project status:

Steps 1-3 are collected and run in the 00pipeline.sh script:
1. 01pvadj.py 
  a. Bonferoni-corrects pvalues
  b. Extracts significant (alpha = 0.05) and non-significant variants (pvalue > 0.9)
  c. Removes unused columns for faster performance downstream
2. 02splitsig.py
  a. Splits data into significant and non-significant subsets
3. 03topsig.py
  a. Gets the most significantly associated variants per gene as positive training set
4. 04getposition.py
  a. Reads a GTF file containing positional information
  b. Finds the position of each variant in relation to its associated gene (intergenic, exonic or intronic)
  c. Stores this information in a new column
5. 05getnegatives.py
  a. Reads user-decided criteria for selection of a negative control (gene ID, position, variant type, etc.)
  b. Finds the most similar non-sig variant for each positive variant in the training set
6. 06getsequences.py
  a. Parses the reference chromosome fasta
  b. Finds the genetic sequence of each variant incl. a +-100 nucleotide window.
  c. Stores the extracted sequence in a new column

Next steps --> find a suitable DNA encoder from literature and use it to encode the extracted sequences for machine learning
