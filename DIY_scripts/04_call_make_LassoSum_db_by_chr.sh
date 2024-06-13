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

# make log directory
mkdir -p step_04_MASHR_LassoSum_db_logs

# Use sbatch --wrap to submit the job for each chromosome
# use the #SBATCH header to specify the resources needed for each job
for i in {1..22}
do
    output_file=step_04_MASHR_LassoSum_db_logs/make_MASHR_LassoSum_db_chr${i}.out
    sbatch --nodes=1 \
        --ntasks-per-node=1 \
        --mem=96G \
        --time=2-12:30:00 \
        --output=${output_file} \
        --wrap="Rscript 04_make_MASHR_db_by_chr.R -f ../pipeline_output/LassoSum_outputs -c AFR-EUR --chr ${i} -g ../pipeline_output/gene_annotation_ALLchr.txt -o ../pipeline_output/MASHR_Lassosum_models"
done

# Rscript 04_make_MASHR_db_by_chr.R \
#     -f ../pipeline_output/LassoSum_outputs \
#     -c AFR-EUR \
#     --chr ${i} \
#     -g ../pipeline_output/gene_annotation_ALLchr.txt \
#     -o ../pipeline_output/MASHR_Lassosum_models




