#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=128G                # Request 64GB of memory
#SBATCH --time=4-12:30:00          # Specify run time

# Load the HPC module environment
source pioneer_HPC_module_env.sh

mkdir -p ../pipeline_output/MatrixeQTL

# start time
date

for pop in AFR EUR HISP
do
  for chr in {1..22}
  do
    Rscript 01_run_MatrixeQTL.R \
    -d ../my_input/processed_input/whole_blood_${pop}_chr${chr}/snp_dosage.txt \
    -e ../my_input/processed_input/whole_blood_${pop}_chr${chr}/gene_expression.txt \
    -g ../my_input/processed_input/whole_blood_${pop}_chr${chr}/gene_annotation.txt \
    -t ${pop}_chr${chr}_cis \
    -o ../pipeline_output/MatrixeQTL/ \
    -w 1000000
  done
done

# end time
date

# time passed
secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs