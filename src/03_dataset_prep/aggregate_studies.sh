#!/bin/bash

# Define the output file name
OUTPUT_FILE="iu_dataset.tsv.gz"

# Find all matching files and store them in an array
mapfile -t FILES < <(find . -type f -name "*iu_dataset.tsv.gz")

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No *iu_dataset.tsv.gz files found!"
    exit 1
fi

echo "Found ${#FILES[@]} files. Starting aggregation..."

# 1. Extract the header from the first file and write it to the output
zcat "${FILES[0]}" | head -n 1 | gzip > "$OUTPUT_FILE"

# 2. Loop through ALL files and append their contents (excluding the header)
for file in "${FILES[@]}"; do
    echo "Processing: $file"
    
    # zcat: read the file
    # tail -n +2: skip the first line (header), print everything from line 2 to the end
    # gzip >>: compress the stream and append it to the output file
    zcat "$file" | tail -n +2 | gzip >> "$OUTPUT_FILE"
done

echo "Aggregation complete! Final file saved as $OUTPUT_FILE"
