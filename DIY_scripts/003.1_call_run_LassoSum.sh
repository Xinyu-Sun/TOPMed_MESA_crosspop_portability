#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=128G                # Request 64GB of memory
#SBATCH --time=4-12:30:00          # Specify run time

# Load the HPC module environment
source pioneer_HPC_module_env.sh

genomebuild='hg38'

for chr in {1..22} # 22
do

  # Use `sbatch --wrap` to run Rscript in the same sbatch script
  sbatch --nodes=1 \
          --ntasks-per-node=2 \
          --mem=128G \
          --time=3-12:30:00 \
          --output=step_003_lassosum_chr${chr}.out \
    --wrap="Rscript 003.1_run_LassoSum.R --chr ${chr}"
done