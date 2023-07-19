#!/bin/bash

# Set the path to the directory containing the files
directory_path="/path/to/your/files"

for file in "$directory_path"/*; do
  # Extract the file name without the path
  file_name=$(basename "$file")

  # Extract the sample name (e.g., 21-01876) from the file name
  sample_name=$(echo "$file_name" | grep -oE '^[0-9]+-[0-9]+')

  # Extract the read type (e.g., R1 or R2) from the file name without the .fastq.gz suffix
  read_type=$(echo "$file_name" | grep -oE '_R[12]' | tail -c 2)

  # Create the new file name using the standard illumina format
  new_file_name="${sample_name}_S1_L001_R${read_type}_001.fastq.gz"

  # Rename the file
  mv "$file" "$directory_path/$new_file_name"

  # Optionally, you can print the changes to the console
  echo "Renamed: $file --> $new_file_name"
done
