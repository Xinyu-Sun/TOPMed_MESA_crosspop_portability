#!/bin/bash

# Usage: ./submit_job.sh arg1 arg2

# Define an array of ancestries
declare -a ancestries=("AA" "NHW" "HISP")

mkdir -p logs

# Loop over each ancestry
for ancestry in "${ancestries[@]}"; do
    # Loop over each chromosome (1 to 22)
    for chromosome in {1..22}; do
        # Set the output log name based on current ancestry and chromosome
        output_log="logs/${ancestry}_${chromosome}_output_%j.log"

        # Submit the job with the current ancestry and chromosome
        sbatch --output="$output_log" slurm.call.process_input.sh $ancestry $chromosome
    done
done

echo "All jobs have been submitted."