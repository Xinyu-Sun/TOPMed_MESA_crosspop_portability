#!/bin/bash

#SBATCH --nodes=1                # Request one node
#SBATCH --ntasks-per-node=2     # Request 16 tasks/cores per node
#SBATCH --mem=64G                # Request 64GB of memory
#SBATCH --time=2-12:30:00          # Specify run time

source ~/.bashrc

module load Miniconda3

conda activate dapg_input
# conda info --env
# conda list 

# start time
date

# $1 is the ancestry group (AA, NHW, HISP)
# $2 is the chromosome number (1..22)
# example usage: bash run_preprocess.sh {AA,NHW,HISP} {1..22}
bash call.process_data_for_input.sh $1 $2 > ${1}.log

# end time
date

# time passed
secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs

conda deactivate