#!/bin/bash

BASE_DIR=~/krispmer/new_datasets/k-mer-counts

for dir in "$BASE_DIR"/*; do
    if [ -d "$dir" ]; then
        input_file="$dir/mer_counts_PAM.txt"
        output_file="$dir/lexo_mer_counts_PAM_sorted.txt"

        if [ -f "$input_file" ]; then
            if [ -f "$output_file" ]; then
                echo "Skipping $dir — sorted file already exists."
            else
                echo "Sorting $input_file by k-mer sequence ..."
                sort -k1,1 "$input_file" > "$output_file"
                echo "  -> Saved to $output_file"
            fi
        else
            echo "No mer_counts_PAM.txt found in $dir"
        fi
    fi
done

