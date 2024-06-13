#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=128G                # Request 64GB of memory
#SBATCH --time=4-12:30:00          # Specify run time

# Load the HPC module environment
source pioneer_HPC_module_env.sh

# This script is for Lassosum + Mashr outputs

# mkae directory for logs

# testing
# pop="EUR"
# gene="ENSG00000005812"

# True parameter values
pop=$1
gene=$2

Rscript 005_make_MASHR_covariances_by_gene.R \
  -d ../my_input/processed_input/whole_blood_ \
  -g ../pipeline_output/gene_annotation_ALLchr.txt \
  -t ${pop}_MASHR \
  --gene ${gene} \
  -m ../pipeline_output/MASHR_Lassosum_models/${pop}_MASHR_weights.txt.gz \
  -o ../pipeline_output/MASHR_Lassosum_models \
  -w 1000000



