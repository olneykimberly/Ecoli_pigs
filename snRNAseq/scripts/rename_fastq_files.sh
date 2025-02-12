#!/bin/bash

# Default directory (can be overridden by an argument)
file_dir="../fastq"

# If an argument is provided, use it as the directory
if [[ -n $1 ]]; then
  file_dir="$1"
fi

# Read the mapping file (adjust path to where `rename_map.txt` is stored)
rename_map_file="rename_map.txt"

# Read the mapping file line by line
while read -r pattern replacement; do
  # Find and rename files matching the pattern
  for file in "$file_dir"/*"$pattern"*; do
    if [[ -f $file ]]; then # Check if it's a file
      # Remove the part of the filename from the second underscore up to the last occurrence of `_S`
      new_name=$(echo "$file" | sed -E 's/(JDFSeq_[0-9]+)_[^_]+_S([0-9]+_L[0-9]+_R[12]_001\.fastq\.gz)$/\1_S\2/')
      # Replace the sample name using the mapping file
      new_name="${new_name//$pattern/$replacement}"
      echo "Renaming: $file -> $new_name"
      mv "$file" "$new_name"
    fi
  done
done < "$rename_map_file"
