#!/bin/bash

# Ensure we are in the right starting directory, or exit if it fails
cd results || { echo "Error: Could not enter 'results' directory."; exit 1; }

echo "Starting cleanup of intermediate files..."
echo "Keeping ONLY files ending in '_dataset.tsv.gz'"

# The find command logic:
# 1. Search in the current directory and all subdirectories (.)
# 2. Only look for files (-type f)
# 3. EXCLUDE files ending in _dataset.tsv.gz (! -name "*_dataset.tsv.gz")
# 4. Delete the matches it finds (-delete)

find . -type f ! -name "*_dataset.tsv.gz" -print -delete

echo "Cleanup complete! All intermediate files have been removed."
