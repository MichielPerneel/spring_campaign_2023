#!/bin/bash

# Define the raw data folder and the new data folder
RAW_DIR="/data/gent/vo/001/gvo00125/vsc43619/spring_campaign_2023/data/raw/RNA"
NEW_DIR="/data/gent/vo/001/gvo00125/vsc43619/HN00231237/RawData"

# Loop through files in the raw directory
for file in "$RAW_DIR"/*.fastq.gz; do
    # Extract the base filename
    basename=$(basename "$file")
    
    # Construct the full path to the corresponding file in the new folder
    new_file="$NEW_DIR/$basename"
    
    # Print details for troubleshooting
    echo "Processing file: $file"
    echo "Looking for: $new_file"
    
    # Check if the corresponding file exists in the new directory
    if [ -f "$new_file" ]; then
        # Create a symbolic link in the raw directory
        echo "Found! Creating symlink: $file -> $new_file"
        ln -sf "$new_file" "$file"
    else
        echo "File not found in new directory: $new_file"
    fi
    
    echo "---------------------------------------------"
done
