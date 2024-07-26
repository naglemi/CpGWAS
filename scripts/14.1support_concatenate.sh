#!/bin/bash

# Define the output file path
output_file="/expanse/lustre/projects/jhu152/naglemi/mwas/pheno/14.1-OUT_combined_methylation_a2.csv"

# Clear the output file to ensure it's empty before we start appending data
> "$output_file"

# List of directories containing the CSV files
directories=(
    "/expanse/lustre/projects/jhu152/naglemi/mwas/pheno/caud/out"
    "/expanse/lustre/projects/jhu152/naglemi/mwas/pheno/hippo/out"
    "/expanse/lustre/projects/jhu152/naglemi/mwas/pheno/dlpfc/out"
)

# Loop through each directory and concatenate all CSV files in it to the output file
for dir in "${directories[@]}"; do
    for file in "$dir"/*_raw-DNAm-only.csv; do
        # Check if file exists to avoid trying to cat non-existent glob expansions
        if [ -f "$file" ]; then
            cat "$file" >> "$output_file"
        else
            echo "No files found in $dir"
        fi
    done
done

echo "Concatenation complete, output in $output_file"

