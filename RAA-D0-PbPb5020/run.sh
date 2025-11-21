#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --job-name=xing
#SBATCH --gres=cpuonly:1

#SBATCH --error=./cpu_%j_%N.err.log
#SBATCH --output=./cpu_%j_%N.out.log

#SBATCH --partition=short

./driver.sh
