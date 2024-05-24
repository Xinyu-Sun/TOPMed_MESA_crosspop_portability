#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=128G                # Request 64GB of memory
#SBATCH --time=4-12:30:00          # Specify run time

# Load the HPC module environment
source pioneer_HPC_module_env.sh

# make output directory for MASHR outputs
mkdir -p ../pipeline_output/MASHR_outputs

# gzip all files under MASHR_inputs if there are any txt file
for file in ../pipeline_output/MASHR_inputs/*.txt
do
  gzip $file
done

# start time
date

for chr in {1..22}
do
  # Rscript 03_run_MASHR.R \
  # -i ../pipeline_output/MASHR_inputs \
  # -g ../my_input/processed_input/whole_blood_EUR_chr${chr}/gene_annotation.txt \
  # -o ../pipeline_output/MASHR_outputs

  # Use `sbatch --wrap` to run Rscript in the same sbatch script
  sbatch --nodes=1 \
          --ntasks-per-node=2 \
          --mem=128G \
          --time=4-12:30:00 \
          --output=step_03_MASHR_chr${chr}.out \
    --wrap="Rscript 03_run_MASHR.R -i ../pipeline_output/MASHR_inputs -g ../my_input/processed_input/whole_blood_EUR_chr${chr}/gene_annotation.txt -o ../pipeline_output/MASHR_outputs"
done

# end time
date

# calculate runtime in H:M:S
secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
