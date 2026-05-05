#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main(file_path):

    print(f"Loading dataset: {file_path}")
    
    # Dropped median_tpm as requested
    columns_to_load = [
        'chromosome', 'variant_location', 
        'type', 'sig_ge', 'sig_iu'
    ]
    
    df = pd.read_csv(file_path, sep='\t', compression='gzip', usecols=columns_to_load)
    
    # ---------------------------------------------------------
    # Prepare Data Categories
    # ---------------------------------------------------------
    # Create a unified 'Effect Type' column
    conditions = [
        (df['sig_ge'] == 1) & (df['sig_iu'] == 1),
        (df['sig_ge'] == 1) & (df['sig_iu'] == 0),
        (df['sig_ge'] == 0) & (df['sig_iu'] == 1)
    ]
    choices = ['Both (GE & IU)', 'Expression (GE) Only', 'Isoform (IU) Only']
    df['Effect_Type'] = np.select(conditions, choices, default='Neither')

    # Format Chromosomes to ensure they start with 'chr' for clean plotting
    df['chr_name'] = df['chromosome'].astype(str)
    df['chr_name'] = df['chr_name'].apply(lambda x: x if x.lower().startswith('chr') else 'chr' + x)
    
    # Define custom chromosome order
    unique_chrs = df['chr_name'].unique().tolist()
    spec_chrs = ['chr1', 'chr19', 'chr8', 'chr2', 'chr5', 'chr16']
    
    # Extract the remaining chromosomes and sort them numerically
    def chr_sort_key(c):
        val = c.replace('chr', '').replace('CHR', '')
        return int(val) if val.isdigit() else 999

    rest_chrs = [c for c in unique_chrs if c not in spec_chrs]
    rest_chrs.sort(key=chr_sort_key)
    
    chr_order = spec_chrs + rest_chrs

    # ---------------------------------------------------------
    # Presentation-Ready Visualizations
    # ---------------------------------------------------------
    print("\nGenerating presentation-ready plots...")
    
    sns.set_theme(context="talk", style="whitegrid")
    
    # Your custom color scheme
    custom_palette = {
        'Both (GE & IU)': '#02882F',        # Saturated Green
        'Expression (GE) Only': '#025A88',  # Mellow Blue
        'Isoform (IU) Only': '#008A74',     # Mellow Turquoise
        'Neither': '#051985'                # Saturated Blue
    }
    effect_order = ['Both (GE & IU)', 'Expression (GE) Only', 'Isoform (IU) Only', 'Neither']

    # --- Plot 1: Global Variant Effects (Donut Chart) ---
    plt.figure(figsize=(8, 8))
    effect_counts = df['Effect_Type'].value_counts().reindex(effect_order)
    
    plt.pie(effect_counts, labels=effect_counts.index, autopct='%1.1f%%', 
            startangle=90, colors=[custom_palette[k] for k in effect_counts.index],
            wedgeprops=dict(width=0.4, edgecolor='black'),  # Black edges
            textprops={'color': 'white', 'weight': 'bold'}) # White text for labels/percentages
            
    plt.title('Global Distribution of Variant Effects', weight='bold', pad=20, color='white')
    plt.tight_layout()
    plt.savefig('results/output/dataset_prep/fig1_global_effects.png', dpi=300, transparent=True)
    plt.close()

    # --- Plot 2: Variant Location vs Effect (100% Stacked Bar) ---
    plt.figure(figsize=(10, 6))
    loc_crosstab = pd.crosstab(df['variant_location'], df['Effect_Type'], normalize='index') * 100
    loc_crosstab = loc_crosstab[effect_order]
    
    loc_crosstab.plot(kind='bar', stacked=True, figsize=(10, 6), 
                      color=[custom_palette[c] for c in loc_crosstab.columns])
    plt.title('Variant Effect Proportions by Genomic Location', weight='bold', pad=15)
    plt.ylabel('Percentage (%)')
    plt.xlabel('Genomic Location')
    plt.xticks(rotation=0)
    plt.legend(title='Effect Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('results/output/dataset_prep/fig2_location_vs_effect.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()

    # --- Plot 3: Chromosome Distribution (Simple Barplot) ---
    plt.figure(figsize=(14, 6))
    chr_counts = df['chr_name'].value_counts().reindex(chr_order)
    
    sns.barplot(x=chr_counts.index, y=chr_counts.values, color='#025A88') # Uses your mellow blue
    plt.title('Total Variants per Chromosome', weight='bold', pad=15)
    plt.ylabel('Variant Count')
    plt.xlabel('Chromosome')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('results/output/dataset_prep/fig3_chromosome_distribution.png', dpi=300, transparent=True)
    plt.close()

    # --- Plot 4: Effect Distribution Across Chromosomes (100% Stacked Bar) ---
    plt.figure(figsize=(16, 6))
    chr_effect_cross = pd.crosstab(df['chr_name'], df['Effect_Type'], normalize='index') * 100
    chr_effect_cross = chr_effect_cross.reindex(chr_order)[effect_order]
    
    chr_effect_cross.plot(kind='bar', stacked=True, figsize=(16, 6), width=0.8,
                          color=[custom_palette[c] for c in chr_effect_cross.columns])
    plt.title('Proportional Variant Effects Across Chromosomes', weight='bold', pad=15)
    plt.ylabel('Percentage (%)')
    plt.xlabel('Chromosome')
    plt.xticks(rotation=45)
    plt.legend(title='Effect Type', bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('results/output/dataset_prep/fig4_chromosome_vs_effect.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()

    # --- Plot 5: Variant Type (SNP vs INDEL) Total Distribution & Breakdown ---
    plt.figure(figsize=(9, 6))
    # Note: Not using 'normalize' here so it plots raw absolute counts
    type_crosstab = pd.crosstab(df['type'], df['Effect_Type'])
    type_crosstab = type_crosstab[effect_order]
    
    type_crosstab.plot(kind='bar', stacked=True, figsize=(9, 6), 
                       color=[custom_palette[c] for c in type_crosstab.columns])
    plt.title('Total Count of SNPs and INDELs by Effect', weight='bold', pad=15)
    plt.ylabel('Total Variant Count')
    plt.xlabel('Variant Type')
    plt.xticks(rotation=0)
    # Format y-axis with commas for readability (e.g., 100,000)
    plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    plt.legend(title='Effect Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('results/output/dataset_prep/fig5_type_absolute_counts.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()

    print("Successfully saved 5 presentation-ready plots to the current directory.")

if __name__ == "__main__":
    file_path = "results/output/dataset_prep/final_full_dataset.tsv.gz"
    main(file_path)
