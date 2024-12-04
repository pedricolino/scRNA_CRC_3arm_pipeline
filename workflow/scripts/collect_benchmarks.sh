#!/bin/bash

benchmark_dir='benchmarks'
target_file='results/stats/benchmarks.tsv'
mkdir -p 'results/stats'

# Remove the target file if it exists
if [ -f "$target_file" ]; then
    rm "$target_file"
fi

# List to store all dataframes
dfs=()

# Traverse through all subdirectories and files
files=$(find "$benchmark_dir" -type f)

for file_path in $files; do
    # Skip the target file
    if [ "$file_path" == "$target_file" ]; then
        continue
    fi
    # Read the file into a DataFrame
    df=$(awk 'BEGIN {FS="\t"; OFS="\t"} NR==1 {print $0, "file_path", "rule"} NR>1 {gsub(/^benchmarks\//, "", FILENAME); print $0, FILENAME, gensub(/([^\/\.]+).*/, "\\1", "g", FILENAME)}' "$file_path")
    if [ ${#dfs[@]} -eq 0 ]; then
        # First file, include header
        df=$(awk 'BEGIN {FS="\t"; OFS="\t"} NR==1 {print $0, "file_path", "rule"} NR>1 {gsub(/^benchmarks\//, "", FILENAME); print $0, FILENAME, gensub(/([^\/\.]+).*/, "\\1", "g", FILENAME)}' "$file_path")
    else
        # Subsequent files, exclude header
        df=$(awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {gsub(/^benchmarks\//, "", FILENAME); print $0, FILENAME, gensub(/([^\/\.]+).*/, "\\1", "g", FILENAME)}' "$file_path")
    fi
    # Append the DataFrame to the list
    dfs+=("$df")
done

# Concatenate all dataframes
printf "%s\n" "${dfs[@]}" > "$target_file"