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
    mpiprocs=( 1 2 4 6 8 10 12 14 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 )
    resolutions=( 125 500 2000 4000 )
    folder="out/datacluster"
    mkdir -p $folder    
else  # if on local machine
    folder="out/datalocal"
    mkdir -p $folder    
    mpiprocs=( 1 2 3 4)
    resolutions=( 125 500 2000 )
fi

iterations=10


for resolution in "${resolutions[@]}"
do  
    for procs in "${mpiprocs[@]}"
    do  
        mpirun -n $procs ./out/build/jacobiMPI $resolution $iterations |& tee "./${folder}/jacobiMPI_${resolution}_${iterations}_n_${procs}.log"
        mpirun -n $procs ./out/build/jacobiMPI_float $resolution $iterations |& tee "./${folder}/jacobiMPI_float_${resolution}_${iterations}_n_${procs}.log"
    done
    ./out/build/jacobiSERIAL $resolution $iterations |& tee "./${folder}/jacobiSERIAL_${resolution}_${iterations}_n_${procs}.log"
done

