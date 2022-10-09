#!/bin/bash
#SBATCH --mail-type=NONE 
#SBATCH --time=48:00:00
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1


module load slurm_setup
module load r/4.1.2-gcc11-mkl

srun Rscript $HOME/mkdir/Intervention/Parallel/03/setup.R
