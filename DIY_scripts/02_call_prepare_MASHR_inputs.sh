#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=128G                # Request 64GB of memory
#SBATCH --time=3-12:30:00          # Specify run time

# Load the HPC module environment
source pioneer_HPC_module_env.sh

mkdir -p ../pipeline_output/MASHR_inputs

for chr in {1..22}
do
  Rscript 02_prepare_MASHR_inputs.R \
  -l ../pipeline_output/02_list_of_input/chr${chr}_input_files.txt \
  -g ../my_input/processed_input/whole_blood_EUR_chr${chr}/gene_annotation.txt \
  -c ${chr} \
  -o ../pipeline_output/MASHR_inputs
done
