#!/bin/bash

###########################################################################################
#job name
#SBATCH --job-name=image
#job stdout file
#SBATCH --output=%a_cont.out
#job stderr file
#SBATCH --error=%a_cont.err
#maximum job time in D-HH:MM
#SBATCH --time=30-00:00
#number of parallel processes (tasks) you are requesting - maps to MPI processes
#SBATCH --ntasks=40
#memory per process in MB
#SBATCH --mem-per-cpu=4000
#partition
#SBATCH --partition=standard
#array job
#SBATCH --array=1-117
#############################################################################################

#module load R/4.0.3

echo "Command = 'Rscript dip_silver.R $SLURM_ARRAY_TASK_ID'"

 Rscript dip_silver.R $SLURM_ARRAY_TASK_ID
