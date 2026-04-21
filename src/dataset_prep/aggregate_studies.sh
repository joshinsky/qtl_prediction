#!/bin/bash

# Define the output file name
OUTPUT_FILE="iu_dataset.tsv.gz"

# Find all matching files and store them in an array
cd results/output/sequence_extraction
mapfile -t FILES < <(find . -type f -name "*iu_dataset.tsv.gz")

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No *iu_dataset.tsv.gz files found!"
    exit 1
fi

echo "Found ${#FILES[@]} files. Starting aggregation..."

# Extract the header from the first file and write it to the output
zcat "${FILES[0]}" | head -n 1 | gzip > "$OUTPUT_FILE"

# Loop through all files and append their contents (excluding the header)
for file in "${FILES[@]}"; do
    echo "Processing: $file"
    zcat "$file" | tail -n +2 | gzip >> "$OUTPUT_FILE"
done

echo "Aggregation complete! Final file saved as $OUTPUT_FILE"

mv "$OUTPUT_FILE" "../dataset_prep"
