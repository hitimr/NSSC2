#!/bin/bash
# Number of cores
#SBATCH -c 1
# Runtime of this jobs is less then 10 minutes
#            (hh:mm:ss)
#SBATCH --time=00:10:00
# Clear the environment
module purge > /dev/null 2>&1
# Set OMP_NUM_THREADS to the same value as -c
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# You can start several programs with one script file/submission
./jacobi 256 40
./jacobi 512 40
./jacobi 1024 40
./jacobi 2048 40
#!/bin/bash
#SBATCH --job-name=exercise1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive 
#SBATCH --time=01:00:00


# Clear the environment
module purge > /dev/null 2>&1

# Set OMP_NUM_THREADS to the same value as -c
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# You can start several programs with one script file/submission
./jacobi 256 40
./jacobi 512 40
./jacobi 1024 40
./jacobi 2048 40
