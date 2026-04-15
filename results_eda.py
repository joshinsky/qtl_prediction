import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gzip

def analyze_variants(filepath):
    print(f"Starting analysis of {filepath}...")
    
    # Initialize dictionaries and counters to store our aggregated metrics
    total_variants = 0
    location_counts = {}
    type_counts = {}
    sig_counts = {'significant': 0, 'non_significant': 0}
    chrom_counts = {}

    # Read the gzipped file in chunks of 1 million rows to prevent OOM errors
    chunk_size = 1000000
    
    # Use pandas to stream the data
    with pd.read_csv(filepath, sep='\t', compression='gzip', chunksize=chunk_size, 
                     usecols=['chromosome', 'type', 'variant_location', 'significant', 'non_significant']) as reader:
        
        for i, chunk in enumerate(reader):
            # Update total count
            total_variants += len(chunk)
            
            # 1. Variant Locations (intronic, exonic, intergenic, etc.)
            locs = chunk['variant_location'].value_counts().to_dict()
            for k, v in locs.items():
                location_counts[k] = location_counts.get(k, 0) + v
                
            # 2. Variant Types (SNP vs INDEL)
            types = chunk['type'].value_counts().to_dict()
            for k, v in types.items():
                type_counts[k] = type_counts.get(k, 0) + v
                
            # 3. Significance
            sig_counts['significant'] += chunk['significant'].sum()
            sig_counts['non_significant'] += chunk['non_significant'].sum()
                
            # 4. Chromosome Distribution
            chroms = chunk['chromosome'].value_counts().to_dict()
            for k, v in chroms.items():
                chrom_counts[k] = chrom_counts.get(k, 0) + v
                
            print(f"Processed chunk {i+1} ({(i+1)*1}M rows)")

    # Print Text Summaries
    print("\n" + "="*40)
    print(f"TOTAL VARIANTS ANALYZED: {total_variants:,}")
    print("="*40)
    
    print("\n--- Variant Locations ---")
    for k, v in sorted(location_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"{k}: {v:,} ({(v/total_variants)*100:.1f}%)")
        
    print("\n--- Variant Types ---")
    for k, v in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"{k}: {v:,} ({(v/total_variants)*100:.1f}%)")
        
    print("\n--- Significance ---")
    for k, v in sig_counts.items():
        print(f"{k}: {v:,} ({(v/total_variants)*100:.1f}%)")

    # ----- Plotting -----
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Locations
    pd.Series(location_counts).sort_values().plot(kind='barh', ax=axes[0,0], color='skyblue')
    axes[0,0].set_title('Variant Locations')
    axes[0,0].set_xlabel('Count')
    
    # Plot 2: Types
    pd.Series(type_counts).plot(kind='pie', ax=axes[0,1], autopct='%1.1f%%', colors=['lightcoral', 'lightgreen'])
    axes[0,1].set_title('Variant Types')
    axes[0,1].set_ylabel('')
    
    # Plot 3: Significance
    pd.Series(sig_counts).plot(kind='bar', ax=axes[1,0], color=['royalblue', 'salmon'])
    axes[1,0].set_title('Significance')
    axes[1,0].tick_params(axis='x', rotation=0)
    
    # Plot 4: Chromosome Distribution
    # Sort chromosomes naturally (chr1, chr2 ... chr22, chrX, chrY)
    chrom_series = pd.Series(chrom_counts)
    chrom_series.index = chrom_series.index.str.replace('chr', '')
    
    def sort_chroms(x):
        try: return int(x)
        except: return ord(x[0]) if isinstance(x, str) and len(x)>0 else 999
        
    sorted_idx = sorted(chrom_series.index, key=sort_chroms)
    chrom_series = chrom_series.reindex(sorted_idx)
    
    chrom_series.plot(kind='bar', ax=axes[1,1], color='mediumpurple')
    axes[1,1].set_title('Chromosome Distribution')
    axes[1,1].set_xlabel('Chromosome')
    axes[1,1].set_ylabel('Variant Count')
    
    plt.tight_layout()
    plt.savefig('variant_eda_dashboard.png', dpi=300)
    print("\nPlots saved to 'variant_eda_dashboard.png'")

if __name__ == "__main__":
    # Point this to your actual file
    analyze_variants("ge_dataset.tsv.gz")
