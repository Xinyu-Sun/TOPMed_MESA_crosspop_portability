#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=128G                # Request 64GB of memory
#SBATCH --time=4-12:30:00          # Specify run time

# Load the HPC module environment
source pioneer_HPC_module_env.sh

# This script is for Lassosum + Mashr outputs

# mkae directory for logs
mkdir -p step_05_MASHR_covariances_logs
for pop in AFR EUR
do
  for chr in {1..22}
  do
    # Rscript DIY_scripts/05_make_MASHR_covariances.R \
    # -d ../my_input/processed_input/whole_blood_${pop}_chr${chr}/snp_dosage.txt \
    # -g ../my_input/processed_input/whole_blood_${pop}_chr${chr}/gene_annotation.txt \
    # -t ${pop}_MASHR \
    # -c ${chr} \
    # -m ../pipeline_output/MASHR_Lassosum_models/${pop}_MASHR_weights.txt.gz \
    # -o ../pipeline_output/MASHR_Lassosum_models \
    # -w 1000000

    # Use `sbatch --wrap` to run Rscript in the same sbatch script
    sbatch --nodes=1 \
            --ntasks-per-node=2 \
            --mem=128G \
            --time=2-12:30:00 \
            --output=step_05_MASHR_covariances_logs/step_05_MASHR_covariances_${pop}_chr${chr}.out \
      --wrap="Rscript 05_make_MASHR_covariances.R -d ../my_input/processed_input/whole_blood_${pop}_chr${chr}/snp_dosage.txt -g ../my_input/processed_input/whole_blood_${pop}_chr${chr}/gene_annotation.txt -t ${pop}_MASHR -c ${chr} -m ../pipeline_output/MASHR_Lassosum_models/${pop}_MASHR_weights.txt.gz -o /scratch/users/xxs410 -w 1000000"
  done
done

# # Gather all job IDs of currently running and pending jobs of the user
# dependency=$(squeue -u $USER -t pending,running -ho %i | paste -sd "," -)

# # Condition to check if there are any dependent jobs
# if [ -n "$dependency" ]; then
#     echo "Setting job dependency to wait for jobs: $dependency"
#     dependency_option="--dependency=afterany:$dependency"
# else
#     echo "No current jobs found. Submitting without dependency."
#     dependency_option=""
# fi
