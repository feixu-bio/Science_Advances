#!/bin/bash

#SBATCH --job-name=radial                   #This is the name of your job
#SBATCH --cpus-per-task=4                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=6:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=output.o%j     #These are the STDOUT and STDERR files
#SBATCH --error=error.e%j

#You selected an array of jobs from 1 to 16 with 16 simultaneous jobs
#SBATCH --array=1-40%10
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=fei.xu@unibas.ch        #You will be notified via email when your task ends or fails




cd /scicore/home/mangos/cokiny55/LADs/early/exp2/

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.

module load MATLAB
matlab -nodisplay -r "radial_distribution_slurm(${SLURM_ARRAY_TASK_ID})"
