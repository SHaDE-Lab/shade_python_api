#!/bin/bash

#SBATCH -N 1                              # number of nodes
#SBATCH -c 8                              # number of cores
#SBATCH -t 7-00:00:00                     # time in d-hh:mm:ss
#SBATCH --mem=25G
#SBATCH -p general                        # partition
#SBATCH -q public                         # QOS
#SBATCH -o slurm.%j.out                   # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL                   # Send an e-mail when a job starts, stops, or fails
#SBATCH --mail-user="wkhan17@asu.edu"
#SBATCH --export=NONE                     # Purge the job-submitting shell environment


#Activate our environment
source activate cool_routes_venv

#Run the software/python script
python /home/wkhan17/Thermal-Comfort-Routing/shade_python_api/server/routing.py
