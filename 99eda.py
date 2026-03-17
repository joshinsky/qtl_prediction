#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_and_inspect(file_path):
    """Loads the TSV file and prints basic information."""
    print(f"Loading data from {file_path}...")
    df = pd.read_csv(file_path, sep='\t', compression='gzip')
    
    print("\n--- Basic Information ---")
    print(df.info())
    
    print("\n--- Missing Values ---")
    missing = df.isnull().sum()
    print(missing[missing > 0] if missing.any() else "No missing values found.")
    
    print("\n--- Numerical Summary ---")
    print(df.describe())
    
    return df

def plot_categorical_distributions(df):
    """Plots distributions of categorical genomic features."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Variant Type (SNP vs INDEL)
    sns.countplot(data=df, x='type', ax=axes[0], palette='Set2')
    axes[0].set_title('Distribution of Variant Types')
    
    # Variant Location
    sns.countplot(data=df, y='variant_location', ax=axes[1], palette='Set3', 
                  order=df['variant_location'].value_counts().index)
    axes[1].set_title('Variant Locations')
    
    # Chromosome Distribution
    # Sorting chromosomes numerically if they follow 'chr1', 'chr2' format
    chrom_order = sorted(df['chromosome'].unique(), 
                         key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else 99)
    sns.countplot(data=df, x='chromosome', ax=axes[2], palette='viridis', order=chrom_order)
    axes[2].set_title('Variants per Chromosome')
    axes[2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.show()

def plot_numerical_distributions(df):
    """Plots distributions of p-values and expression levels."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # -log10 Adjusted P-value Distribution
    # Adding a small constant to avoid log(0)
    df['neg_log10_padj'] = -np.log10(df['p_adj'] + 1e-300)
    sns.histplot(df['neg_log10_padj'], bins=50, kde=True, ax=axes[0], color='crimson')
    axes[0].set_title('Distribution of -log10(Adjusted P-value)')
    axes[0].set_xlabel('-log10(p_adj)')
    
    # Median TPM Distribution (Log scaled for better visualization)
    df['log1p_tpm'] = np.log1p(df['median_tpm'])
    sns.histplot(df['log1p_tpm'], bins=50, kde=True, ax=axes[1], color='dodgerblue')
    axes[1].set_title('Distribution of log1p(Median TPM)')
    axes[1].set_xlabel('log1p(TPM)')
    
    plt.tight_layout()
    plt.show()

def plot_significance_vs_expression(df):
    """Scatter plot of Expression (TPM) vs Significance (-log10 P-value)."""
    plt.figure(figsize=(8, 6))
    
    # Separate significant and non-significant for coloring
    sns.scatterplot(
        data=df, 
        x='log1p_tpm', 
        y='neg_log10_padj', 
        hue='significant',
        palette={0: 'grey', 1: 'red'},
        alpha=0.6,
        edgecolor=None
    )
    
    # Add a horizontal line for standard significance threshold (e.g., p_adj < 0.05)
    plt.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p_adj = 0.05')
    
    plt.title('Gene Expression vs Variant Significance')
    plt.xlabel('log1p(Median TPM)')
    plt.ylabel('-log10(Adjusted P-value)')
    plt.legend(title='Significant (1=Yes)')
    plt.tight_layout()
    plt.show()

def analyze_sequences(df):
    """Basic analysis of the variant_window sequence lengths."""
    df['window_length'] = df['variant_window'].apply(len)
    
    plt.figure(figsize=(6, 4))
    sns.histplot(df['window_length'], bins=30, color='purple')
    plt.title('Distribution of Variant Window Sequence Lengths')
    plt.xlabel('Sequence Length (bp)')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.show()
    
    print("\n--- Sequence Window Lengths ---")
    print(df['window_length'].describe())

def main(file_path):
    df = load_and_inspect(file_path)
    plot_categorical_distributions(df)
    plot_numerical_distributions(df)
    plot_significance_vs_expression(df)
    analyze_sequences(df)
    return df

# Run the EDA pipeline
if __name__ == "__main__":
    tsv_file = sys.argv[1]
    eda_df = main(tsv_file)
