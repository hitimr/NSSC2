#!/bin/bash
#SBATCH -J jacobiMPI
#SBATCH -N 2
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive 
#SBATCH --time=01:00:00
if command -v sinfo  2>/dev/null # if on cluster
then
    module load mpi/openmpi-x86_64
    module load pmi/pmix-x86_64
    mpiprocs=( 1 2 5 10 20 40 80 )
    folder="out/datacluster"
    mkdir -p $folder    
else  # if on local machine
    folder="out/datalocal"
    mkdir -p $folder    
    mpiprocs=( 1 2 )
fi

iterations=10
resolutions=( 125 500 2000 4000 )

for resolution in "${resolutions[@]}"
do  
    for procs in "${mpiprocs[@]}"
    do  
        mpirun -n $procs ./out/build/jacobiMPI $resolution $iterations |& tee "./${folder}/jacobiMPI_${resolution}_${iterations}_n_${procs}.log"
        mpirun -n $procs ./out/build/jacobiMPI_float $resolution $iterations |& tee "./${folder}/jacobiMPI_float_${resolution}_${iterations}_n_${procs}.log"
    done
done

