#!/bin/bash
#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=128G                # Request 64GB of memory
#SBATCH --time=4-12:30:00          # Specify run time
#SBATCH --output=make_MASHR_LassoSum_db.out

# Load the HPC module environment
source pioneer_HPC_module_env.sh

# make output directory
mkdir -p ../pipeline_output/MASHR_Lassosum_models

# concate all gene_annotation.txt from {1..22} chromosomes under my_input/processed_input/whole_blood_EUR_chr{1..22}/gene_annotation.txt
# remove the excessive header lines
# save the concated gene_annotation.txt to pipeline_output/gene_annotation_ALLchr.txt
# below is the **WRONG code**
# cat ../my_input/processed_input/whole_blood_EUR_chr{1..22}/gene_annotation.txt > ../pipeline_output/gene_annotation_ALLchr.txt

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


# start time
date

Rscript 04_make_MASHR_db.R \
-f ../pipeline_output/LassoSum_outputs \
-c AFR-EUR \
-g ../pipeline_output/gene_annotation_ALLchr.txt \
-o ../pipeline_output/MASHR_Lassosum_models

# end time
date

# time passed
secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
