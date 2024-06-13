#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 2 tasks/cores per node
#SBATCH --mem=128G                # Request 128GB of memory
#SBATCH --time=4-12:30:00          # Specify run time

# Load the HPC module environment
source pioneer_HPC_module_env.sh

# This script is for Lassosum + Mashr outputs

# Make directory for logs
mkdir -p step_05_MASHR_covariances_logs

# Function to submit a job with dependencies
submit_job() {
  local dependency=$1
  local pop=$2
  local chr=$3

  if [ -z "$dependency" ]; then
    # Submit job without dependency
    sbatch --nodes=1 \
           --ntasks-per-node=2 \
           --mem=128G \
           --time=2-12:30:00 \
           --output=step_05_MASHR_covariances_logs/step_05_MASHR_covariances_${pop}_chr${chr}.out \
           --wrap="Rscript 05_make_MASHR_covariances.R -d ../my_input/processed_input/whole_blood_${pop}_chr${chr}/snp_dosage.txt -g ../my_input/processed_input/whole_blood_${pop}_chr${chr}/gene_annotation.txt -t ${pop}_MASHR -c ${chr} -m ../pipeline_output/MASHR_Lassosum_models/${pop}_MASHR_weights.txt.gz -o ../pipeline_output/MASHR_Lassosum_models -w 1000000"
  else
    # Submit job with dependency
    sbatch --dependency=afterany:$dependency \
           --nodes=1 \
           --ntasks-per-node=2 \
           --mem=128G \
           --time=2-12:30:00 \
           --output=step_05_MASHR_covariances_logs/step_05_MASHR_covariances_${pop}_chr${chr}.out \
           --wrap="Rscript 05_make_MASHR_covariances.R -d ../my_input/processed_input/whole_blood_${pop}_chr${chr}/snp_dosage.txt -g ../my_input/processed_input/whole_blood_${pop}_chr${chr}/gene_annotation.txt -t ${pop}_MASHR -c ${chr} -m ../pipeline_output/MASHR_Lassosum_models/${pop}_MASHR_weights.txt.gz -o ../pipeline_output/MASHR_Lassosum_models -w 1000000"
  fi
}

# Initialize an array to store the job IDs
job_ids=()

# Submit the first 11 jobs without dependencies
counter=0
for pop in AFR EUR; do
  for chr in {1..22}; do
    if [ $counter -lt 11 ]; then
      job_id=$(submit_job "" $pop $chr | awk '{print $NF}')
      job_ids+=($job_id)
      counter=$((counter + 1))
    fi
  done
done

# Submit the remaining jobs with dependencies
for pop in AFR EUR; do
  for chr in {1..22}; do
    if [ $counter -ge 11 ]; then
      dependency_index=$(( (counter - 11) % 11 ))
      dependency_job=${job_ids[$dependency_index]}
      job_id=$(submit_job $dependency_job $pop $chr | awk '{print $NF}')
      job_ids[$dependency_index]=$job_id
      counter=$((counter + 1))
    fi
  done
done
