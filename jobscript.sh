#use this to log in to Vera : ssh patriku@vera2.c3se.chalmers.se 

#!/urs/bin/env bash
#SBATCH -A SNIC2021-7-120
#SBATCH -t 12:00:00
#SBATCH --job-name= Na6 cluster program    #This is just a name so we can recognize it on the queue list 
#SBATCH -N 1 # Use 1 node
#SBATCH -n 1 # Use only 1 core on that node


model purge #this clears the environment from previous modules. 
module load python3
module load intel/2019a GPAW ASE

#run the code:
gpaw-python ga.py
#./ga.py ? 
#srun ./ga.py ? 

# skicka in detta i terminalen för att submitta jobbet : sbatch jobscript.sh
# squeue JOB ID to see information 


########################################################################################################################

#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2 # Project
#SBATCH -J GPAW # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 1 # Use only 1 core on that node
#SBATCH -t 10:00:00 # Maximum time
#SBATCH -o std.out # stdout goes to this file
#SBATCH -e err.out # stderr goes to this file

# Unload all modules, to make sure no incompatibilities arise
module purge

# Load desired modules
module load GPAW/22.8.0-foss-2022a

# Run program
srun gpaw python ga.py

##sbatch jobscript.sh #this we put in the terminal to actually submit the job. 
