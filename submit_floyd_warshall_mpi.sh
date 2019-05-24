#!/bin/bash
#
#SBATCH --job-name=floyd-warshall-50-nodes
#SBATCH --nodes=50
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4GB
#SBATCH --mail-type=END
#SBATCH --mail-user=asa566@nyu.edu
#SBATCH --output=results/mpi_50_nodes_results.out

module purge
module load gcc/6.3.0
module load openmpi/intel/3.1.3

cd /scratch/asa566/hpc/HPC-final-project
mpirun -np 50 ./floyd-warshall-mpi
