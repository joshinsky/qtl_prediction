#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    file_path = "results/output/dataset_prep/final_full_dataset.tsv.gz"
    print(f"Loading dataset: {file_path}")
    
    columns_to_load = [
        'molecular_trait_id', 'chromosome', 'variant_location', 
        'type', 'sig_ge', 'sig_iu', 'median_tpm'
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

    print("\n" + "="*50)
    print(f"TOTAL DATAPOINTS: {len(df):,}")
    print("="*50)
    print(df['Effect_Type'].value_counts(normalize=True) * 100)

    # ---------------------------------------------------------
    # Presentation-Ready Visualizations
    # ---------------------------------------------------------
    print("\nGenerating presentation-ready plots with custom color scheme...")
    
    # Set seaborn context to 'talk' for larger presentation fonts
    sns.set_theme(context="talk", style="whitegrid")
    
    # Apply your custom presentation color scheme
    custom_palette = {
        'Both (GE & IU)': '#02882F',        # Saturated Green
        'Expression (GE) Only': '#025A88',  # Mellow Blue
        'Isoform (IU) Only': '#008A74',     # Mellow Turquoise
        'Neither': '#051985'                # Saturated Blue
    }

    # Ordering for legends and axes
    effect_order = ['Both (GE & IU)', 'Expression (GE) Only', 'Isoform (IU) Only', 'Neither']

    # --- Plot 1: Global Variant Effects (Donut Chart) ---
    plt.figure(figsize=(8, 8))
    effect_counts = df['Effect_Type'].value_counts().reindex(effect_order)
    
    plt.pie(effect_counts, labels=effect_counts.index, autopct='%1.1f%%', 
            startangle=90, colors=[custom_palette[k] for k in effect_counts.index],
            wedgeprops=dict(width=0.4, edgecolor='w'))
    plt.title('Global Distribution of Variant Effects', weight='bold', pad=20)
    plt.tight_layout()
    plt.savefig('results/figures/presentation/fig1_global_effects.png', dpi=300, transparent=True)
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
    plt.savefig('results/figures/presentation/fig2_location_vs_effect.png', dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 3: Variant Type (SNP vs INDEL) vs Effect ---
    plt.figure(figsize=(8, 6))
    type_crosstab = pd.crosstab(df['type'], df['Effect_Type'], normalize='index') * 100
    type_crosstab = type_crosstab[effect_order]
    
    type_crosstab.plot(kind='bar', stacked=True, figsize=(8, 6), 
                       color=[custom_palette[c] for c in type_crosstab.columns])
    plt.title('Variant Effect Proportions by Variant Type', weight='bold', pad=15)
    plt.ylabel('Percentage (%)')
    plt.xlabel('Variant Type')
    plt.xticks(rotation=0)
    plt.legend(title='Effect Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('results/figures/presentation/fig3_type_vs_effect.png', dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 4: Expression Level Bias (Boxplot with Log Scale) ---
    plt.figure(figsize=(10, 6))
    df_expr = df[df['median_tpm'] > 0]
    
    sns.boxplot(data=df_expr, x='Effect_Type', y='median_tpm', 
                order=effect_order, palette=custom_palette)
    
    plt.yscale('log')
    plt.title('Gene Expression Levels (TPM) across Effect Types', weight='bold', pad=15)
    plt.ylabel('Median TPM (Log Scale)')
    plt.xlabel('Variant Effect Type')
    plt.xticks(rotation=15)
    plt.tight_layout()
    plt.savefig('results/figures/presentation/fig4_tpm_bias.png', dpi=300)
    plt.close()

    print("Successfully saved 4 presentation-ready plots to the current directory:")
    print(" - results/figures/presentation/fig1_global_effects.png")
    print(" - results/figures/presentation/fig2_location_vs_effect.png")
    print(" - results/figures/presentation/fig3_type_vs_effect.png")
    print(" - results/figures/presentation/fig4_tpm_bias.png")

if __name__ == "__main__":
    main()
