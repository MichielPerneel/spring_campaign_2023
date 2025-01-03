#!/bin/bash

# Run from the spring_campaign_2023 directory on the HPC

# Define directories for the two runs
RUN1_DIR="../HN00231237/RawData"
RUN2_DIR="../EN00004671"

# Define the output directory for combined files
OUTPUT_DIR="data/raw/RNA_seqruns_combined"
mkdir -p "${OUTPUT_DIR}"

# Loop over 130_* and 51_* samples from the first run
for sample in "${RUN1_DIR}"/130_*_R1.fastq.gz "${RUN1_DIR}"/51_*_R1.fastq.gz; do
    # Extract the sample name without the _R1.fastq.gz suffix
    sample_name=$(basename "${sample}" _R1.fastq.gz)

    # Define the corresponding R1 and R2 files for both runs
    run1_r1="${RUN1_DIR}/${sample_name}_R1.fastq.gz"
    run1_r2="${RUN1_DIR}/${sample_name}_R2.fastq.gz"
    run2_r1="${RUN2_DIR}/${sample_name}_1.fastq.gz"
    run2_r2="${RUN2_DIR}/${sample_name}_2.fastq.gz"

    # Define output file paths
    combined_r1="${OUTPUT_DIR}/${sample_name}_combined_R1.fastq.gz"
    combined_r2="${OUTPUT_DIR}/${sample_name}_combined_R2.fastq.gz"

    echo "Processing ${sample_name}..."

    # Check if R1 files exist and concatenate them
    if [[ -f "${run1_r1}" && -f "${run2_r1}" ]]; then
        echo "Combining ${sample_name} R1 files from both runs..."
        cat "${run1_r1}" "${run2_r1}" > "${combined_r1}"
    elif [[ -f "${run1_r1}" ]]; then
        echo "Second run R1 missing for ${sample_name}, using only first run R1."
        cp "${run1_r1}" "${combined_r1}"
    else
        echo "Warning: First run R1 missing for ${sample_name}, skipping."
        continue
    fi

    # Check if R2 files exist and concatenate them
    if [[ -f "${run1_r2}" && -f "${run2_r2}" ]]; then
        echo "Combining ${sample_name} R2 files from both runs..."
        cat "${run1_r2}" "${run2_r2}" > "${combined_r2}"
    elif [[ -f "${run1_r2}" ]]; then
        echo "Second run R2 missing for ${sample_name}, using only first run R2."
        cp "${run1_r2}" "${combined_r2}"
    else
        echo "Warning: First run R2 missing for ${sample_name}, skipping."
        continue
    fi

    echo "Combined files saved to ${combined_r1} and ${combined_r2}"
done

echo "Combining complete!"
