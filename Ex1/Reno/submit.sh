#!/bin/bash
#SBATCH --job-name=exercise1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --exclusive 
#SBATCH --time=01:00:00

if command -v sinfo  2>/dev/null # if on cluster
then
    module load mpi/openmpi-x86_64
    module load pmi/pmix-x86_64
    mpiprocs=( 1 2 3 4 5 6 7 8 9 10 
11 12 13 14 15 16 17 18 19 20 
21 22 23 24 25 26 27 28 29 30 
31 32 33 34 35 36 37 38 39 40 
41 42 43 44 45 46 47 48 49 50 
51 52 53 54 55 56 57 58 59 60 
61 62 63 64 65 66 67 68 69 70
 71 72 73 74 75 76 77 78 79 80 )

    folder="datacluster"
    mkdir -p $folder    
else  # if on local machine
    folder="datalocal"
    mkdir -p $folder    
    mpiprocs=( 1 2 )
fi

iterations=10
resolutions=( 125 50 1000 2000 )

for resolution in "${resolutions[@]}"
do  
    for procs in "${mpiprocs[@]}"
    do  
        mpirun -n $procs ./jacobiMPI $resolution $iterations |& tee "./${folder}/jacobiMPI_${resolution}_${iterations}_n_${procs}.log"
    done
done