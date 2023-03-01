#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2 # Project
#SBATCH -J GPAW # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 30 # Use only 1 core on that node
#SBATCH -t 23:00:00 # Maximum time
#SBATCH -o std.out # stdout goes to this file
#SBATCH -e err.out # stderr goes to this file
# Unload all modules, to make sure no incompatibilities arise
module purge

# Load desired modules
module load GPAW/22.8.0-foss-2022a

# Run program
srun gpaw python task6.py

##sbatch jobscript.sh #this we put in the terminal to actually submit the job.
