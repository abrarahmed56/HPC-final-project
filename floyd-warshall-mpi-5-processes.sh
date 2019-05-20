#!/bin/bash
#
##SBATCH --nodes=5
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=1:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=floyd-warshall-5-processes
#SBATCH --mail-type=END
#SBATCH --mail-user=asa566@nyu.edu
#SBATCH --output=5_processes_results_updated_code.out

module purge
module load gcc/6.3.0
module load openmpi/intel/3.1.3

cd /scratch/asa566/hpc/HPC-final-project
./floyd-warshall-mpi
