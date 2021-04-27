#!/bin/bash

#SBATCH --job-name=int-ring
#SBATCH --nodes=2 
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB 
#SBATCH --time=00:30:00 
#SBATCH --verbose

module purge 
module load  amber/openmpi/intel/20.06 

mpirun -np 25 --oversubscribe ./int_ring 100
