#!/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def apply_dark_slide_aesthetics(ax):
    """Helper function to format axes for dark presentation slides."""
    ax.set_facecolor('none') # Transparent axes background
    ax.tick_params(colors='white', width=1.5)
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')
    # Make the outer bounding box of the graph bold and white
    for spine in ax.spines.values():
        spine.set_color('white')
        spine.set_linewidth(1.5)

def style_legend(ax):
    """Helper function to format the legend for dark slides."""
    leg = ax.legend(title='Effect Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.setp(leg.get_texts(), color='white')
    plt.setp(leg.get_title(), color='white', weight='bold')
    leg.get_frame().set_alpha(0.0) # Transparent legend box
    leg.get_frame().set_edgecolor('white')

def main(in_file_path, out_file_path):

    os.makedirs(out_file_path, exist_ok=True)
    print(f"Loading dataset: {in_file_path}")
    
    columns_to_load = [
        'chromosome', 'variant_location', 
        'type', 'sig_ge', 'sig_iu'
    ]
    
    df = pd.read_csv(in_file_path, sep='\t', compression='gzip', usecols=columns_to_load)
    
    # ---------------------------------------------------------
    # Prepare Data Categories
    # ---------------------------------------------------------
    conditions = [
        (df['sig_ge'] == 1) & (df['sig_iu'] == 1),
        (df['sig_ge'] == 1) & (df['sig_iu'] == 0),
        (df['sig_ge'] == 0) & (df['sig_iu'] == 1)
    ]
    choices = ['Both (GE & IU)', 'Expression (GE) Only', 'Isoform (IU) Only']
    df['Effect_Type'] = np.select(conditions, choices, default='Neither')

    df['chr_name'] = df['chromosome'].astype(str)
    df['chr_name'] = df['chr_name'].apply(lambda x: x if x.lower().startswith('chr') else 'chr' + x)
    
    unique_chrs = df['chr_name'].unique().tolist()
    spec_chrs = ['chr1', 'chr19', 'chr8', 'chr2', 'chr5', 'chr16']
    
    def chr_sort_key(c):
        val = c.replace('chr', '').replace('CHR', '')
        return int(val) if val.isdigit() else 999

    rest_chrs = [c for c in unique_chrs if c not in spec_chrs]
    rest_chrs.sort(key=chr_sort_key)
    chr_order = spec_chrs + rest_chrs

    effect_order = ['Both (GE & IU)', 'Expression (GE) Only', 'Isoform (IU) Only', 'Neither']
    
    custom_palette = {
        'Both (GE & IU)': '#02882F',        # Saturated Green
        'Expression (GE) Only': '#025A88',  # Mellow Blue
        'Isoform (IU) Only': '#008A74',     # Mellow Turquoise
        'Neither': '#051985'                # Saturated Blue
    }

    # ---------------------------------------------------------
    # Terminal Printouts
    # ---------------------------------------------------------
    print("\n" + "="*50)
    print("VARIANT LOCATION STATISTICS")
    print("="*50)
    print("Absolute Counts:")
    print(df['variant_location'].value_counts())
    
    loc_crosstab = pd.crosstab(df['variant_location'], df['Effect_Type'], normalize='index') * 100
    loc_crosstab = loc_crosstab[effect_order]
    print("\nPercentages of Effect Types per Location:")
    print(loc_crosstab.round(2))

    print("\n" + "="*50)
    print("VARIANT TYPE (SNP/INDEL) STATISTICS")
    print("="*50)
    type_crosstab = pd.crosstab(df['type'], df['Effect_Type'])
    type_crosstab = type_crosstab[effect_order]
    type_crosstab_pct = type_crosstab.div(type_crosstab.sum(axis=1), axis=0) * 100
    print("Percentages of Effect Types per Variant Type:")
    print(type_crosstab_pct.round(2))


    # ---------------------------------------------------------
    # Presentation-Ready Visualizations
    # ---------------------------------------------------------
    print("\nGenerating presentation-ready plots...")
    
    # Transparent global figure settings
    sns.set_theme(context="talk", style="whitegrid", rc={"axes.facecolor": (0, 0, 0, 0), "figure.facecolor": (0, 0, 0, 0), 'grid.color': '#444444'})

    # --- Plot 1: Global Variant Effects (Donut Chart) ---
    plt.figure(figsize=(8, 8))
    effect_counts = df['Effect_Type'].value_counts().reindex(effect_order)
    
    plt.pie(effect_counts, labels=effect_counts.index, autopct='%1.1f%%', 
            startangle=90, colors=[custom_palette[k] for k in effect_counts.index],
            wedgeprops=dict(width=0.4, edgecolor='black', linewidth=1.5),  # Bold black edges
            textprops={'color': 'white', 'weight': 'bold'}) # White text
            
    plt.title('Global Distribution of Variant Effects', weight='bold', pad=20, color='white')
    plt.tight_layout()
    plt.savefig(f'{out_file_path}/fig1_global_effects.png', dpi=300, transparent=True)
    plt.close()

    # --- Plot 2: Variant Location vs Effect (100% Stacked Bar) ---
    plt.figure(figsize=(10, 6))
    ax = loc_crosstab.plot(kind='bar', stacked=True, figsize=(10, 6), 
                           color=[custom_palette[c] for c in loc_crosstab.columns],
                           edgecolor='black', linewidth=1.5) # Bold black edges
                           
    # Add percentage labels to bars
    pct_values = loc_crosstab.values.T.flatten()
    for i, p in enumerate(ax.patches):
        width, height = p.get_width(), p.get_height()
        if height > 1.5:  # Only label segments larger than 1.5% to avoid overlapping text
            ax.text(p.get_x() + width/2, p.get_y() + height/2, f'{pct_values[i]:.1f}%',
                    ha='center', va='center', color='white', weight='bold', fontsize=12)

    apply_dark_slide_aesthetics(ax)
    plt.title('Variant Effect Proportions by Genomic Location', weight='bold', pad=15)
    plt.ylabel('Percentage (%)')
    plt.xlabel('Genomic Location')
    plt.xticks(rotation=0)
    style_legend(ax)
    plt.tight_layout()
    plt.savefig(f'{out_file_path}/fig2_location_vs_effect.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()

    # --- Plot 3: Chromosome Distribution (Simple Barplot) ---
    plt.figure(figsize=(14, 6))
    chr_counts = df['chr_name'].value_counts().reindex(chr_order)
    
    ax = sns.barplot(x=chr_counts.index, y=chr_counts.values, color='#025A88', edgecolor='black', linewidth=1.5)
    apply_dark_slide_aesthetics(ax)
    plt.title('Total Variants per Chromosome', weight='bold', pad=15)
    plt.ylabel('Variant Count')
    plt.xlabel('Chromosome')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{out_file_path}/fig3_chromosome_distribution.png', dpi=300, transparent=True)
    plt.close()

    # --- Plot 4: Effect Distribution Across Chromosomes (100% Stacked Bar) ---
    plt.figure(figsize=(16, 6))
    chr_effect_cross = pd.crosstab(df['chr_name'], df['Effect_Type'], normalize='index') * 100
    chr_effect_cross = chr_effect_cross.reindex(chr_order)[effect_order]
    
    ax = chr_effect_cross.plot(kind='bar', stacked=True, figsize=(16, 6), width=0.8,
                               color=[custom_palette[c] for c in chr_effect_cross.columns],
                               edgecolor='black', linewidth=1.5)
    apply_dark_slide_aesthetics(ax)
    plt.title('Proportional Variant Effects Across Chromosomes', weight='bold', pad=15)
    plt.ylabel('Percentage (%)')
    plt.xlabel('Chromosome')
    plt.xticks(rotation=45)
    style_legend(ax)
    plt.tight_layout()
    plt.savefig(f'{out_file_path}/fig4_chromosome_vs_effect.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()

    # --- Plot 5: Variant Type (SNP vs INDEL) Total Distribution & Breakdown ---
    plt.figure(figsize=(9, 6))
    ax = type_crosstab.plot(kind='bar', stacked=True, figsize=(9, 6), 
                            color=[custom_palette[c] for c in type_crosstab.columns],
                            edgecolor='black', linewidth=1.5)
                            
    # Add percentage labels to the absolute count bars
    pct_values_type = type_crosstab_pct.values.T.flatten()
    for i, p in enumerate(ax.patches):
        width, height = p.get_width(), p.get_height()
        if height > (type_crosstab.values.sum() * 0.005):  # Only label significant vertical chunks
            ax.text(p.get_x() + width/2, p.get_y() + height/2, f'{pct_values_type[i]:.1f}%',
                    ha='center', va='center', color='white', weight='bold', fontsize=12)

    apply_dark_slide_aesthetics(ax)
    plt.title('Total Count of SNPs and INDELs by Effect', weight='bold', pad=15)
    plt.ylabel('Total Variant Count')
    plt.xlabel('Variant Type')
    plt.xticks(rotation=0)
    plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    style_legend(ax)
    plt.tight_layout()
    plt.savefig(f'{out_file_path}/fig5_type_absolute_counts.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()

    print("Successfully saved 5 presentation-ready plots to the current directory.")

if __name__ == "__main__":
    in_file_path = "results/output/dataset_prep/final_full_dataset.tsv.gz"
    out_file_path = "results/figures/presentation"
    main(in_file_path, out_file_path)
